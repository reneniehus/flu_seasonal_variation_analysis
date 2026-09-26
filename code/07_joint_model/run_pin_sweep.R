# run_pin_sweep.R -- what does the R0 pin do? Evidence for the owner's decision (MODEL.md, "The R0 pin").
#
# R0 and S0 are exactly one number to the data: (k R0, S0/k, I0/k, k c) reproduces every expected count.
# So the pin can matter through three things only, and this sweep measures each:
#   1. the ceiling S0 <= 1, which binds at a low pin (fit, cells near 1);
#   2. the logit-additive form of the season and country effects, which a rescaling of S0 does not
#      preserve (fit, and whether rankings and season patterns move);
#   3. the ABSOLUTE scale: at a fixed R0 x S0 the infections scale by 1/k and the reporting fraction by k,
#      so the pin sets the attack rate -- the one place where outside evidence (serology) can inform it.
# The S0 prior centre moves with the pin (0.75 x 1.5 / R0), so the prior sits on the same R0 x S0 at
# every pin and cannot decide between them. About 12 minutes on four cores.
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_compare.R"); jm_load_cpp()
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj
options(width = 200)
out = "output/joint_model/pin"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
cands = c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")
cores = max(1L, parallel::detectCores() - 1L)
pins = c(1.5, 1.7, 2.0, 2.5, 3.0)
fits = lapply(pins, function(r0){
  f = file.path(out, sprintf("fit_R0_%.1f.rds", r0))
  if (file.exists(f)) return(readRDS(f))
  s = modifyList(jm_settings(), list(R0_fixed = r0, pr_S0_mean = qlogis(0.75 * 1.5 / r0)))
  d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  fit = jm_fit(d, theta0 = jm_theta0(d, s), cores = cores, verbose = FALSE)
  saveRDS(fit, f); fit
})
names(fits) = sprintf("%.1f", pins)
ref = fits[["1.5"]]; p_ref = jm_unpack(ref$theta, ref$d)
s0c = function(p) vapply(p$country, `[[`, 0, "S0")
rows = lapply(names(fits), function(k){
  fit = fits[[k]]; d = fit$d; p = jm_unpack(fit$theta, d); f = jm_fitted_cpp(fit$theta, d)
  s0 = unlist(lapply(p$country, `[[`, "S0_season"))
  # population-weighted attack rate per country-season, and by age
  att = t(vapply(seq_len(d$n_cs), function(i){ N = d$N[[d$cs_country[i] + 1L]]
    c(all = sum(f$attack[i, ] * N) / sum(N), f$attack[i, ]) }, numeric(4)))
  data.frame(R0 = d$R0_fixed, loglik = fit$loglik, d_loglik = fit$loglik - ref$loglik, converged = fit$converged,
             S0_median = median(s0), S0_max = max(s0), cells_above_0.95 = sum(s0 > 0.95),
             Reff_median = median(d$R0_fixed * s0),
             S0_rank_vs_1.5 = cor(s0c(p), s0c(p_ref), method = "spearman"),
             season_rank_vs_1.5 = cor(p$x, p_ref$x, method = "spearman"),
             visibility_r_vs_1.5 = cor(p$delta, p_ref$delta),
             attack_all = median(att[, "all"]), attack_young = median(att[, 2]), attack_adult = median(att[, 3]),
             attack_elderly = median(att[, 4]),
             reporting_adult_median = median(vapply(p$country, `[[`, 0, "c")), row.names = NULL)
})
tab = do.call(rbind, rows)
cat("=== the R0 pin: fit, ceiling, what moves, and the absolute scale it sets ===\n")
print(tab, row.names = FALSE, digits = 3)
cmp = do.call(rbind, lapply(names(fits)[-1], function(k)
  jm_compare_waves(ref, fits[[k]], sprintf("R0 %s vs 1.5", k))))
cat("\n=== wave-level comparison against the pin of 1.5 ===\n")
print(cmp[, c("contrast", "favour_B", "favour_A", "sign_p", "wilcoxon_p", "total", "total_deflated",
              "ci_wave_lo", "ci_wave_hi", "ci_country_lo", "ci_country_hi")], row.names = FALSE, digits = 3)
saveRDS(list(table = tab, compare = cmp), file.path(out, "pin_sweep.rds"))
cat("DONE\n")
