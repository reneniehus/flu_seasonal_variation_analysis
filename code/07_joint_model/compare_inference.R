# compare_inference.R -- inference for the model comparison that respects WHAT IS INDEPENDENT.
#
# Weekly observations inside one wave are autocorrelated: the deterministic mean does not follow the
# week-to-week wander (the noise budget measures it), so residuals of adjacent weeks are correlated. A
# raw log-likelihood difference sums over those weeks as if each were an independent observation, and
# so does AIC. Both therefore OVERSTATE the evidence, by roughly the effective-sample-size factor.
#
# The unit of evidence here is the WAVE (country-season): 85 of them. Everything below is built on
# per-wave paired differences in log-likelihood, and uses only between-wave independence:
#   - sign test and Wilcoxon signed-rank on the per-wave differences;
#   - a cluster bootstrap of the TOTAL difference, resampling whole waves -- and, because the shared
#     season and country parameters couple waves, also whole countries and whole seasons, which is
#     more conservative;
#   - the residual lag-1 autocorrelation inside waves, so the raw nats can be shown next to an
#     independent-equivalent scale (a heuristic deflation, reported as such -- the bootstrap is the
#     primary inference).
suppressMessages({library(dplyr)})

# Pearson residuals per wave under a fit, and their lag-1 autocorrelation
.jm_resid_acf1 = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d); f = jm_fitted_cpp(fit$theta, d)
  vapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; y = d$y[[i]]; mu = f$mu[[i]]; phi = p$country[[ic]]$phi
    # pool the age groups by summing counts: one series per wave, which is what the dependence is
    # about; within-week cross-age dependence is a separate matter
    yt = rowSums(y); mt = rowSums(mu); ok = is.finite(yt) & rowSums(is.finite(y)) == ncol(y)
    if (sum(ok) < 8) return(NA_real_)
    r = (yt[ok] - mt[ok]) / sqrt(mt[ok] + mt[ok]^2 / phi)
    r = r - mean(r); n = length(r)
    sum(r[-1] * r[-n]) / sum(r^2)
  }, numeric(1))
}

# one contrast: model B minus model A, per wave, with wave-level inference
jm_compare_contrast = function(cells, A, B, label, B_boot = 4000L, seed = 1L){
  a = cells[cells$model == A, ]; b = cells[cells$model == B, ]
  m = merge(a, b, by = c("cell", "country", "season"), suffixes = c("_A", "_B"))
  dll = m$loglik_B - m$loglik_A                      # > 0 favours B
  n = length(dll); n_pos = sum(dll > 0); n_neg = sum(dll < 0)
  sign_p = binom.test(n_pos, n_pos + n_neg, 0.5)$p.value
  wilc_p = suppressWarnings(wilcox.test(dll, mu = 0, exact = FALSE)$p.value)
  set.seed(seed)
  boot = function(groups){                           # resample whole groups, sum their dll
    g = split(dll, groups); k = length(g)
    tot = replicate(B_boot, sum(unlist(g[sample.int(k, k, replace = TRUE)])))
    quantile(tot, c(0.025, 0.975))
  }
  ci_wave = boot(m$cell); ci_country = boot(m$country); ci_season = boot(m$season)
  data.frame(contrast = label, A = A, B = B, n_waves = n, favour_B = n_pos, favour_A = n_neg,
             sign_p = sign_p, wilcoxon_p = wilc_p, total = sum(dll), median_per_wave = median(dll),
             ci_wave_lo = ci_wave[1], ci_wave_hi = ci_wave[2],
             ci_country_lo = ci_country[1], ci_country_hi = ci_country[2],
             ci_season_lo = ci_season[1], ci_season_hi = ci_season[2],
             row.names = NULL, stringsAsFactors = FALSE)
}

jm_compare_inference = function(cmp, fits, B_boot = 4000L){
  cells = cmp$cells
  contrasts = list(
    c("S0/S0", "R0/S0", "season effect: R0 vs S0  (country on S0)"),
    c("S0/S0", "S0/R0", "country effect: R0 vs S0  (season on S0)"),
    c("S0/R0", "R0/R0", "season effect: R0 vs S0  (country on R0)"),
    c("R0/S0", "R0/R0", "country effect: R0 vs S0  (season on R0)"),
    c("S0/S0", "R0/R0", "both on R0 vs both on S0"))
  tab = do.call(rbind, lapply(contrasts, function(k)
    jm_compare_contrast(cells, k[1], k[2], k[3], B_boot = B_boot)))
  # residual autocorrelation inside waves, per model, and the implied effective-sample-size factor
  acf = do.call(rbind, lapply(names(fits), function(m){
    r = .jm_resid_acf1(fits[[m]])
    data.frame(model = m, waves = sum(is.finite(r)), acf1_median = median(r, na.rm = TRUE),
               acf1_q25 = quantile(r, .25, na.rm = TRUE), acf1_q75 = quantile(r, .75, na.rm = TRUE),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  rho = median(acf$acf1_median)
  ess_factor = (1 - rho) / (1 + rho)                 # AR(1) effective-sample-size ratio, a heuristic
  tab$total_deflated = tab$total * ess_factor
  list(contrasts = tab, acf = acf, rho = rho, ess_factor = ess_factor)
}

if (sys.nframe() == 0){
  setwd(here::here())
  suppressMessages(source("code/01_main_supporting/setup.R"))
  source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
  source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
  source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
  source("code/07_joint_model/joint_model.R"); jm_load_cpp()
  options(width = 220)
  cmp = readRDS("output/joint_model/compare/comparison.rds")
  fits = list(); for (m in c("S0/S0", "R0/R0", "R0/S0", "S0/R0"))
    fits[[m]] = readRDS(sprintf("output/joint_model/compare/fit_%s.rds", gsub("[/ =.]", "_", m)))
  inf = jm_compare_inference(cmp, fits)
  cat("================ residual lag-1 autocorrelation INSIDE waves (why weekly nats overstate) ================\n")
  print(inf$acf, row.names = FALSE, digits = 3)
  cat(sprintf("\n  median rho = %.3f  ->  AR(1) effective-sample-size factor (1-rho)/(1+rho) = %.2f\n",
              inf$rho, inf$ess_factor))
  cat("  so a raw difference of N nats is worth roughly N x that factor in independent-equivalent terms\n")
  cat("\n================ the contrasts, WAVE as the unit (85 country-seasons) ================\n")
  t = inf$contrasts
  for (i in seq_len(nrow(t))) with(t[i, ], {
    cat(sprintf("\n%s\n", contrast))
    cat(sprintf("  waves favouring %s: %d ; favouring %s: %d   sign test p = %.3g ; Wilcoxon p = %.3g\n",
                B, favour_B, A, favour_A, sign_p, wilcoxon_p))
    cat(sprintf("  total %+.1f nats (median per wave %+.2f) ; independent-equivalent ~ %+.1f\n",
                total, median_per_wave, total_deflated))
    cat(sprintf("  95%% cluster-bootstrap CI on the total: waves [%+.1f, %+.1f] ; countries [%+.1f, %+.1f] ; seasons [%+.1f, %+.1f]\n",
                ci_wave_lo, ci_wave_hi, ci_country_lo, ci_country_hi, ci_season_lo, ci_season_hi))
  })
  saveRDS(inf, "output/joint_model/compare/inference.rds")
  cat("\nDONE\n")
}
