# comp_model_fit.R -- per-country fit of the compartmental model (ASSUMPTIONS.md sections 8, 9)
#
# Parameters per country (transformed scale for optim):
#   logit_S0            one initial susceptibility, shared across seasons and age groups (E2)
#   log_R0[1..K]        per-season transmissibility (only if R0_free; else fixed at settings$R0_reference)
#   log_c               reporting proportion (fraction of observation-relevant infections that are ILI+)
#   log_b               off-season baseline, ILI+ rate per 100 000 (converted to counts per group)
#   log_phi             observation overdispersion (Var = mu + mu^2/phi)
#   log_q               multiplicative process-noise sd on the infected fractions
# Penalties (priors): logit S0 ~ N(logit .75, 1); log R0_s ~ N(log 1.5, 0.05); log phi ~ N(log 15, .8);
# log q ~ N(log .05, 1); c, b unpenalised. Multi-start BFGS via .fit_multistart() (sir_core.R).

cm_unpack = function(theta, K, R0_free, settings){
  i = 1
  S0 = plogis(theta[i]); i = i + 1
  if (R0_free){ R0 = exp(theta[i:(i + K - 1)]); i = i + K } else R0 = rep(settings$R0_reference, K)
  list(S0 = S0, R0 = R0, c = exp(theta[i]), b = exp(theta[i + 1]), phi = exp(theta[i + 2]), q = exp(theta[i + 3]))
}

cm_logprior = function(theta, K, R0_free, settings){
  i = 1; lp = dnorm(theta[i], settings$prior_logitS0["mean"], settings$prior_logitS0["sd"], log = TRUE); i = i + 1
  if (R0_free){ lp = lp + sum(dnorm(theta[i:(i + K - 1)], log(settings$R0_reference), settings$prior_logR0_sd, log = TRUE)); i = i + K }
  lp + dnorm(theta[i + 2], settings$prior_logphi["mean"], settings$prior_logphi["sd"], log = TRUE) +
       dnorm(theta[i + 3], settings$prior_logq["mean"],   settings$prior_logq["sd"],   log = TRUE)
}

# ---- |-penalised negative log-likelihood over a country's seasons ----
cm_negll = function(theta, cd, f, settings, R0_free = TRUE, return_fit = FALSE){
  K = length(cd$seasons); p = cm_unpack(theta, K, R0_free, settings)
  engine = if (is.null(settings$engine) || !exists("cm_ekf_season_engine", mode = "function")) "R" else settings$engine
  ll = 0; filt = vector("list", K)
  for (s in seq_len(K)){
    e = if (engine == "R") cm_ekf_season(cd$y[[s]], f, p$S0, p$R0[s], p$c, p$b, p$phi, p$q, cd$vax_day, cd$vax_frac[[s]])
        else cm_ekf_season_engine(cd$y[[s]], f, p$S0, p$R0[s], p$c, p$b, p$phi, p$q, cd$vax_day, cd$vax_frac[[s]], engine = engine)
    ll = ll + e$loglik; filt[[s]] = e
  }
  lp = cm_logprior(theta, K, R0_free, settings)
  if (return_fit) return(list(negll = -(ll + lp), loglik = ll, logprior = lp, params = p, filt = filt))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
}

# ---- |-fit one country ----
fit_comp_model = function(cd, settings, R0_free = TRUE, n_starts = settings$n_starts, seed = 1, verbose = TRUE,
                          cores = max(1, min(n_starts, parallel::detectCores() - 1))){
  f = cm_fixed(cd$Cn, cd$N, settings); K = length(cd$seasons)
  # data-driven starts: reporting from the observed peak (peak weekly incidence ~2% of a group), baseline from the floor
  rate_all = unlist(lapply(cd$rates, function(m) as.numeric(m)))
  peak_rate = max(rate_all, na.rm = TRUE); floor_rate = max(as.numeric(quantile(rate_all[rate_all > 0], 0.1, na.rm = TRUE)), 1e-3)
  c0 = min(max(peak_rate / settings$rate_per / 0.02, 1e-4), 1)     # mu = c N C -> c ~ (peak rate/1e5) / peak C
  base = c(qlogis(0.8), if (R0_free) rep(log(settings$R0_reference), K), log(c0), log(floor_rate), log(15), log(0.05))
  names(base) = c("logit_S0", if (R0_free) paste0("log_R0_", cd$seasons), "log_c", "log_b", "log_phi", "log_q")
  jit = c(0.6, if (R0_free) rep(0.03, K), 0.5, 0.5, 0.5, 0.6)
  best = .fit_multistart(base, jit, cm_negll,
                         list(cd = cd, f = f, settings = settings, R0_free = R0_free),
                         n_starts, seed, sprintf("fit_comp_model(%s)", cd$country),
                         cores = cores, maxit = settings$optim_maxit)
  fit = cm_negll(best$par, cd, f, settings, R0_free, return_fit = TRUE)
  # curvature-based standard errors on the transformed scale (Laplace)
  H = tryCatch(optimHess(best$par, cm_negll, cd = cd, f = f, settings = settings, R0_free = R0_free), error = function(e) NULL)
  se = if (!is.null(H)) tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(best$par))) else rep(NA_real_, length(best$par))
  names(se) = names(base)
  # deterministic curves + attack rates at the optimum (the mechanistic mean, no filtering)
  det = lapply(seq_len(K), function(s){
    sim = cm_simulate_season(f, fit$params$S0, fit$params$R0[s], cd$n_weeks[s], cd$vax_day, cd$vax_frac[[s]])
    list(mu = cm_mu(sim$inc, f, fit$params$c, fit$params$b), attack = sim$attack)
  })
  if (verbose) cat(sprintf("[%s] S0=%.3f  R0: %s  c=%.3g b=%.2f phi=%.1f q=%.3f  negll=%.1f conv=%d\n",
                           cd$country, fit$params$S0, paste(sprintf("%.3f", fit$params$R0), collapse=" "),
                           fit$params$c, fit$params$b, fit$params$phi, fit$params$q, best$value, best$convergence))
  list(country = cd$country, seasons = cd$seasons, groups = cd$groups, N = cd$N, R0_free = R0_free,
       theta = best$par, se = se, params = fit$params, negll = best$value, loglik = fit$loglik,
       convergence = best$convergence, settings = settings,
       y = cd$y, rates = cd$rates, mu_filt = lapply(fit$filt, `[[`, "mu_pred"),
       I_filt = lapply(fit$filt, `[[`, "I"), S_filt = lapply(fit$filt, `[[`, "S"),
       mu_det = lapply(det, `[[`, "mu"), attack = do.call(rbind, lapply(det, `[[`, "attack")),
       R_eff = fit$params$R0 * fit$params$S0, vax = cd$vax, contact_source = cd$contact_source)
}

# ---- |-tidy per-season summary of a fit (one row per season) ----
summarise_comp_fit = function(fit){
  K = length(fit$seasons)
  data.frame(country = fit$country, season = fit$seasons, S0 = fit$params$S0, R0 = fit$params$R0,
             R_eff = fit$R_eff, c = fit$params$c, b = fit$params$b, phi = fit$params$phi, q = fit$params$q,
             attack_young = fit$attack[, 1], attack_medium = fit$attack[, 2], attack_elderly = fit$attack[, 3],
             cor = vapply(seq_len(K), function(s){ y = as.numeric(fit$y[[s]]); m = as.numeric(fit$mu_filt[[s]])
               ok = is.finite(y) & is.finite(m); if (sum(ok) > 2) cor(y[ok], m[ok]) else NA_real_ }, numeric(1)),
             vax_cov_65 = fit$vax$coverage, vax_provenance = fit$vax$provenance, stringsAsFactors = FALSE)
}
