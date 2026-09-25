# run_model_comparison.R -- WHERE DOES THE VARIATION LIVE: in the susceptible pool, or in transmissibility?
#
# Four models with IDENTICAL parameter counts, differing only in which quantity carries the season
# effect and which the country effect:
#
#             country effect on S0        country effect on R0
#   season    S0/S0  the working model    S0/R0  hybrid
#   effect    R0/S0  hybrid (close to the R0/R0
#   on R0            previous model)
#
# Whatever carries no effect is pinned: R0 = 1.5, S0 = 0.75. Same data, same likelihood, same count, so
# the models compare by log-likelihood directly (AIC and BIC differ from it by a constant). They are not
# nested, so this is a likelihood comparison, not a test.
#
# What the two mechanisms predict differently: an effect on S0 changes a wave's rise rate AND the pool
# left to infect (so its final size); an effect on R0 changes the rise rate only. So the per-cell
# likelihood differences say, wave by wave, which shape the data prefer.
#
# Reported: total and per-cell log-likelihood, convergence and flat-line status for each, the
# breakdown by season and by country, a sensitivity of the R0/R0 result to the pinned S0, and the
# prior-to-posterior contraction of R0/R0 so the comparison is known to be between maximum-likelihood
# fits rather than prior-driven ones.
#   Rscript code/07_joint_model/run_model_comparison.R
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_recovery.R")
jm_load_cpp()
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj
cores = max(1L, parallel::detectCores() - 1L)
candidates = c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")
dir.create("output/joint_model/compare", showWarnings = FALSE, recursive = TRUE)

# per-cell log-likelihood, computed in R from the fitted means: exact for the negative binomial and
# the only way to see WHICH waves prefer which mechanism
cell_loglik = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d); f = jm_fitted_cpp(fit$theta, d)
  vapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; y = d$y[[i]]; mu = f$mu[[i]]; ok = is.finite(y)
    sum(dnbinom(y[ok], size = p$country[[ic]]$phi, mu = pmax(mu[ok], 1e-10), log = TRUE))
  }, numeric(1))
}

variants = list(
  `S0/S0` = list(season_on = "S0", country_on = "S0"),
  `R0/R0` = list(season_on = "R0", country_on = "R0"),
  `R0/S0` = list(season_on = "R0", country_on = "S0"),
  `S0/R0` = list(season_on = "S0", country_on = "R0"),
  `R0/R0 S0=0.85` = list(season_on = "R0", country_on = "R0", S0_fixed = 0.85)   # sensitivity to the pin
)

fits = list(); rows = list(); cells = list()
for (nm in names(variants)){
  set = modifyList(jm_settings(), variants[[nm]])
  d = jm_build_data(candidates, models_in, demo, set = set, verbose = FALSE)
  cat(sprintf("\n=== %-14s season on %s, country on %s, R0 = %.2f, S0 = %.2f, %d parameters ===\n",
              nm, set$season_on, set$country_on, set$R0_fixed, set$S0_fixed, d$n_par))
  t0 = Sys.time()
  fit = jm_fit(d, cores = cores, verbose = FALSE)
  secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ll = jm_loglik_cpp(fit$theta, d)
  cl = cell_loglik(fit)
  cat(sprintf("  loglik %.2f  negll %.2f  converged %s  flat %d  %.0f s  (sum of cell logliks %.2f)\n",
              ll, fit$negll, fit$converged, fit$n_flat_unresolved, secs, sum(cl)))
  fits[[nm]] = fit
  rows[[nm]] = data.frame(model = nm, season_on = set$season_on, country_on = set$country_on,
                          R0_fixed = set$R0_fixed, S0_fixed = set$S0_fixed, k = d$n_par,
                          loglik = ll, negll = fit$negll, converged = fit$converged,
                          n_flat = fit$n_flat_unresolved, seconds = secs, stringsAsFactors = FALSE)
  cells[[nm]] = data.frame(model = nm, cell = seq_len(d$n_cs),
                           country = d$countries[d$cs_country + 1L], season = d$seasons[d$cs_season + 1L],
                           loglik = cl, stringsAsFactors = FALSE)
  saveRDS(fit, sprintf("output/joint_model/compare/fit_%s.rds", gsub("[/ =.]", "_", nm)))
}
tab = do.call(rbind, rows); rownames(tab) = NULL
tab$AIC = -2 * tab$loglik + 2 * tab$k
tab$dAIC_vs_S0S0 = tab$AIC - tab$AIC[tab$model == "S0/S0"]
cat("\n================ the comparison ================\n")
print(tab[, c("model", "season_on", "country_on", "S0_fixed", "k", "loglik", "AIC", "dAIC_vs_S0S0", "converged", "n_flat")],
      row.names = FALSE, digits = 6)
cat("\n(reference: the previous model, R0_s fitted with a free mean + S0_c, 182 parameters, loglik -50786.68)\n")

# the 2x2 decomposition: main effects of WHERE the season effect lives and where the country effect lives
g = function(m) tab$loglik[tab$model == m]
cat("\n================ the 2x2, in log-likelihood ================\n")
cat(sprintf("  season effect on R0 instead of S0, averaged over the country choice: %+.2f nats\n",
            mean(c(g("R0/S0") - g("S0/S0"), g("R0/R0") - g("S0/R0")))))
cat(sprintf("  country effect on R0 instead of S0, averaged over the season choice: %+.2f nats\n",
            mean(c(g("S0/R0") - g("S0/S0"), g("R0/R0") - g("R0/S0")))))
cat(sprintf("  interaction (does the answer for one depend on the other?):         %+.2f nats\n",
            (g("R0/R0") - g("R0/S0")) - (g("S0/R0") - g("S0/S0"))))

# per-cell: which waves prefer a season effect on R0 over S0, holding the country choice at S0
cc = merge(cells[["S0/S0"]], cells[["R0/S0"]], by = c("cell", "country", "season"), suffixes = c("_S", "_R"))
cc$dll = cc$loglik_R - cc$loglik_S           # > 0 favours the season effect on R0
cat("\n================ per country-season: season effect on R0 minus on S0 (country on S0 in both) ================\n")
cat(sprintf("  cells favouring R0: %d of %d ; total %+.2f nats ; median per cell %+.3f ; range %+.2f to %+.2f\n",
            sum(cc$dll > 0), nrow(cc), sum(cc$dll), median(cc$dll), min(cc$dll), max(cc$dll)))
bys = aggregate(dll ~ season, cc, sum); byc = aggregate(dll ~ country, cc, sum)
cat("  by season:\n"); print(bys[order(bys$season), ], row.names = FALSE, digits = 3)
cat("  by country:\n"); print(byc[order(-byc$dll), ], row.names = FALSE, digits = 3)
cat("\n  the ten cells that most favour R0 and the ten that most favour S0:\n")
print(head(cc[order(-cc$dll), c("country", "season", "dll")], 10), row.names = FALSE, digits = 3)
print(head(cc[order(cc$dll), c("country", "season", "dll")], 10), row.names = FALSE, digits = 3)

saveRDS(list(table = tab, cells = do.call(rbind, cells), pairwise = cc),
        "output/joint_model/compare/comparison.rds")

# contraction of the R0/R0 model, so the comparison is known to be between data-driven optima
cat("\n================ R0/R0: does the data decide, or the priors? ================\n")
id = jm_identifiability(fits[["R0/R0"]])
print(id$family[, c("family", "n", "prior_sd", "sd_post_med", "contraction_med", "contraction_min")],
      row.names = FALSE, digits = 3)
cat(sprintf("near-flat directions %d, penalised Hessian positive definite %s\n", id$n_near_flat, id$post_pd))
saveRDS(id, "output/joint_model/compare/ident_R0R0.rds")
cat("\nDONE\n")
