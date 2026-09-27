# run_kappa_profile.R -- the children's S0 modifier: what the data say about it, and where its prior
# comes from (MODEL.md, appendix L).
#
# kappa raises children's S0 on the logit scale (logit S0_young = logit S0 + kappa). The surveillance
# data barely separate it from children's reporting, so its prior is anchored on cohort evidence:
# centred on the kappa at which the REFITTED model infects children as often, relative to adults, as the
# PHIRST cohort found (about 1.7; Cohen et al. 2021). This runner holds kappa at a grid of values with a
# needle prior, refits everything else, and reports
#   1. the profile of the fit, and a wave-level comparison of each value against none (joint_compare.R);
#   2. the children-to-adult infection ratio at each value, from which the prior centre is read.
# kappa acts on the logit scale, so the ratio a given kappa produces depends on the level of S0 and with
# it on the R0 pin: re-run this after changing the pin and move pr_kappa_mean to the new calibration.
# About 25 minutes on three cores.
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_compare.R"); jm_load_cpp()
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj
options(width = 200)
out = "output/joint_model/kappa"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
o = readRDS("output/joint_model/joint_fit.rds"); cands = o$fit$d$countries
cores = max(1L, parallel::detectCores() - 1L)
PHIRST = 1.7                                  # children-to-adult infection ratio (decisions.md: 1.65-1.80)
grid = c(0, 0.25, 0.5, 1, 1.5, 2, 3)
fits = lapply(grid, function(k){
  lk = if (k == 0) -30 else log(k)            # exp(-30): the modifier switched off
  s = modifyList(jm_settings(), list(pr_kappa_mean = lk, pr_kappa_sd = 0.005))   # held by a needle prior
  d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  f = file.path(out, sprintf("fit_kappa_%g.rds", k))
  if (file.exists(f)){ fit = readRDS(f); if (isTRUE(all.equal(fit$d, d))) return(fit) }   # same data and settings only
  th0 = o$fit$theta; th0[["log_kappa_young"]] = lk
  fit = jm_fit(d, theta0 = th0, cores = cores, verbose = FALSE); saveRDS(fit, f); fit
})
ref = fits[[1]]
tab = do.call(rbind, Map(function(fit, k){
  d = fit$d; p = jm_unpack(fit$theta, d); a = jm_fitted_cpp(fit$theta, d)$attack
  data.frame(kappa = k, d_loglik = fit$loglik - ref$loglik, converged = fit$converged,
             young_adult = median(a[, 1] / a[, 2]), elderly_adult = median(a[, 3] / a[, 2]),
             S0_adult = median(unlist(lapply(p$country, `[[`, "S0_season"))),
             S0_young = median(unlist(lapply(p$country, `[[`, "S0_young_season"))),
             child_reporting = median(2^vapply(p$country, `[[`, 0, "off_young")), sigma_eld = p$sigma_eld)
}, fits, grid))
cat("=== profile over the children's S0 modifier (everything else refitted; ratios are medians over country-seasons) ===\n")
print(tab, row.names = FALSE, digits = 3)
cmp = Map(function(fit, k) jm_compare_waves(ref, fit, sprintf("kappa %g vs none", k)), fits[-1], grid[-1])
cat("\n=== wave-level comparison against no modifier (positive favours the modifier) ===\n")
print(do.call(rbind, cmp)[, c("contrast", "favour_B", "favour_A", "sign_p", "wilcoxon_p", "total",
                             "ci_wave_lo", "ci_wave_hi", "ci_country_lo", "ci_country_hi")], row.names = FALSE, digits = 3)
by_country = sapply(cmp, function(x){ w = attr(x, "per_wave"); tapply(w$dll, w$country, sum) })
colnames(by_country) = sprintf("kappa %g", grid[-1])
cat("\n=== where the cost falls: log-likelihood difference against no modifier, by country ===\n")
print(round(by_country, 1))
centre = approx(tab$young_adult, tab$kappa, xout = PHIRST)$y
cat(sprintf("\nthe refitted model reproduces PHIRST's children-to-adult ratio of %.1f at kappa = %.2f (prior centre in jm_settings: %.2f)\n",
            PHIRST, centre, exp(o$fit$d$pr_kappa_mean)))
saveRDS(list(table = tab, compare = do.call(rbind, cmp), by_country = by_country, centre = centre),
        file.path(out, "profile_kappa.rds"))
cat("DONE\n")
