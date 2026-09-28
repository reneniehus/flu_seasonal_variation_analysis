# joint_compare.R -- comparing two fits of the joint model with the WAVE as the unit of evidence.
#
# Weekly observations inside one wave are autocorrelated: the deterministic mean does not follow the
# week-to-week wander, so residuals of adjacent weeks are correlated (lag-1 ~0.2). A raw difference in
# log-likelihood sums over those weeks as if each were independent, and so does AIC; both overstate the
# evidence. Everything here is built on PER-WAVE paired differences in log-likelihood and uses only
# between-wave independence:
#   - sign test and Wilcoxon signed-rank on the per-wave differences;
#   - a cluster bootstrap of the total difference, resampling whole waves, whole countries and whole
#     seasons (the shared parameters couple waves, so the coarser clusters are the conservative ones);
#   - the residual lag-1 autocorrelation, so the raw nats can be shown next to an independent-
#     equivalent scale (a heuristic deflation, reported as such; the bootstrap is the inference).
# First used for the S0-versus-R0 sensing comparison (2026-09-25, MODEL.md); kept for every model
# comparison since, because that is what the learning layer does.
suppressMessages({library(dplyr)})

# Per-wave log-likelihood under a fit, from R's own dnbinom -- a computation independent of the C++.
jm_wave_loglik = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d); f = jm_fitted_cpp(fit$theta, d)
  data.frame(cell = seq_len(d$n_cs), country = d$countries[d$cs_country + 1L],
             season = d$seasons[d$cs_season + 1L],
             loglik = vapply(seq_len(d$n_cs), function(i){
               ic = d$cs_country[i] + 1L; y = d$y[[i]]; mu = f$mu[[i]]; ok = is.finite(y)
               sum(dnbinom(y[ok], size = p$country[[ic]]$phi, mu = pmax(mu[ok], 1e-10), log = TRUE))
             }, numeric(1)), stringsAsFactors = FALSE)
}

# Lag-1 autocorrelation of the Pearson residuals inside each wave (age groups pooled by summing).
jm_resid_acf1 = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d); f = jm_fitted_cpp(fit$theta, d)
  vapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; y = d$y[[i]]; mu = f$mu[[i]]; phi = p$country[[ic]]$phi
    yt = rowSums(y); mt = rowSums(mu); ok = is.finite(yt) & rowSums(is.finite(y)) == ncol(y)
    if (sum(ok) < 8) return(NA_real_)
    r = (yt[ok] - mt[ok]) / sqrt(mt[ok] + mt[ok]^2 / phi); r = r - mean(r); n = length(r)
    sum(r[-1] * r[-n]) / sum(r^2)
  }, numeric(1))
}

# Model B against model A, both fitted to the SAME data. Positive differences favour B.
jm_compare_waves = function(fit_a, fit_b, label = "B vs A", B_boot = 4000L, seed = 1L){
  if (!identical(fit_a$d$y, fit_b$d$y)) stop("the two fits are not on the same data")
  a = jm_wave_loglik(fit_a); b = jm_wave_loglik(fit_b)
  dll = b$loglik - a$loglik
  n_pos = sum(dll > 0); n_neg = sum(dll < 0)
  set.seed(seed)
  boot = function(groups){ g = split(dll, groups); k = length(g)
    quantile(replicate(B_boot, sum(unlist(g[sample.int(k, k, replace = TRUE)]))), c(0.025, 0.975)) }
  ci_w = boot(a$cell); ci_c = boot(a$country); ci_s = boot(a$season)
  rho = median(c(jm_resid_acf1(fit_a), jm_resid_acf1(fit_b)), na.rm = TRUE)
  ess = (1 - rho) / (1 + rho)                         # AR(1) effective-sample-size ratio, a heuristic
  out = data.frame(contrast = label, n_waves = length(dll), favour_B = n_pos, favour_A = n_neg,
                   sign_p = binom.test(n_pos, n_pos + n_neg, 0.5)$p.value,
                   wilcoxon_p = suppressWarnings(wilcox.test(dll, mu = 0, exact = FALSE)$p.value),
                   total = sum(dll), median_per_wave = median(dll), rho = rho, total_deflated = sum(dll) * ess,
                   ci_wave_lo = ci_w[[1]], ci_wave_hi = ci_w[[2]],
                   ci_country_lo = ci_c[[1]], ci_country_hi = ci_c[[2]],
                   ci_season_lo = ci_s[[1]], ci_season_hi = ci_s[[2]],
                   k_A = length(fit_a$theta), k_B = length(fit_b$theta), row.names = NULL,
                   stringsAsFactors = FALSE)
  attr(out, "per_wave") = data.frame(a[, c("cell", "country", "season")], dll = dll)
  out
}
