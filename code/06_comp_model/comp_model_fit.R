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

# parameter vector layout:
#   logit_S0 | [log_R0 x K] | [log_I0 x K] | log_c | [log2 c-offset young, elderly] | [log c-deviation x K] | log_b | log_phi | log_q
# c is returned as a K x A matrix (season x age) so every consumer indexes c[s, ] regardless of the switches.
cm_unpack = function(theta, K, R0_free, settings){
  i = 1
  S0 = plogis(theta[i]); i = i + 1
  if (R0_free){ R0 = exp(theta[i:(i + K - 1)]); i = i + K } else R0 = rep(settings$R0_reference, K)
  if (isTRUE(settings$I0_by_season)){ I0 = exp(theta[i:(i + K - 1)]); i = i + K } else I0 = rep(settings$I0_fraction, K)
  c_med = exp(theta[i]); i = i + 1
  if (isTRUE(settings$c_by_age)){ off = theta[i:(i + 1)]; i = i + 2; c_age = c_med * 2^c(off[1], 0, off[2]) }   # young, medium (ref), elderly
  else c_age = rep(c_med, 3)
  if (isTRUE(settings$c_by_season)){ dev = exp(theta[i:(i + K - 1)]); i = i + K } else dev = rep(1, K)
  c_mat = outer(unname(dev), unname(c_age))                                                        # K x A
  if (isTRUE(settings$susc_by_age)){ sg = theta[i:(i + 1)]; i = i + 2; sigma = 2^c(sg[1], 0, sg[2]) }   # young, medium (ref), elderly
  else sigma = if (!is.null(settings$susc_fixed)) settings$susc_fixed else c(1, 1, 1)                  # or an externally FIXED profile
  if (isTRUE(settings$b_by_source)){ b = exp(theta[i:(i + 1)]); names(b) = c("RespiCompass", "ERVISS"); i = i + 2 }
  else b = unname(exp(theta[i])); if (!isTRUE(settings$b_by_source)) i = i + 1
  q = if (is.null(settings$q_fixed)) unname(exp(theta[i + 1])) else settings$q_fixed
  list(S0 = unname(S0), R0 = unname(R0), I0 = unname(I0), c = c_mat, c_age = unname(c_age), c_season = unname(dev),
       sigma = unname(sigma), b = b, phi = unname(exp(theta[i])), q = q)
}

# fixed inputs with the age-susceptibility profile applied: rows of Cn scaled by sigma, renormalised
# to spectral radius 1 (so R0_s is still realised exactly; sigma only redistributes infection by age)
cm_fixed_sigma = function(f, sigma){
  if (all(sigma == 1)) return(f)
  Cs = sweep(f$Cn, 1, sigma, "*"); f$Cn = Cs / spectral_radius(Cs); f
}

# baseline of season s: per source when b_by_source, else the shared scalar
cm_b_season = function(p, cd, s) if (length(p$b) > 1) unname(p$b[cd$source_by_season[s]]) else p$b

cm_logprior = function(theta, K, R0_free, settings){
  i = 1; lp = dnorm(theta[i], settings$prior_logitS0["mean"], settings$prior_logitS0["sd"], log = TRUE); i = i + 1
  if (R0_free){ lp = lp + sum(dnorm(theta[i:(i + K - 1)], log(settings$R0_reference), settings$prior_logR0_sd, log = TRUE)); i = i + K }
  if (isTRUE(settings$I0_by_season)){ lp = lp + sum(dnorm(theta[i:(i + K - 1)], log(settings$I0_fraction), settings$prior_logI0_sd, log = TRUE)); i = i + K }
  if (!is.null(settings$prior_logc)) lp = lp + dnorm(theta[i], settings$prior_logc["mean"], settings$prior_logc["sd"], log = TRUE)
  i = i + 1
  if (isTRUE(settings$c_by_age)){ lp = lp + sum(dnorm(theta[i:(i + 1)], 0, settings$prior_logc_age_sd, log = TRUE)); i = i + 2 }
  if (isTRUE(settings$c_by_season)){ lp = lp + sum(dnorm(theta[i:(i + K - 1)], 0, settings$prior_logc_season_sd, log = TRUE)); i = i + K }
  if (isTRUE(settings$susc_by_age)){ lp = lp + sum(dnorm(theta[i:(i + 1)], 0, settings$prior_log2susc_sd, log = TRUE)); i = i + 2 }
  i = i + if (isTRUE(settings$b_by_source)) 2 else 1                     # log_b: unpenalised
  lp = lp + dnorm(theta[i], settings$prior_logphi["mean"], settings$prior_logphi["sd"], log = TRUE)
  if (is.null(settings$q_fixed)) lp = lp + dnorm(theta[i + 1], settings$prior_logq["mean"], settings$prior_logq["sd"], log = TRUE)
  lp
}

# ---- |-deterministic-model likelihood of one season (stage 1: no filter) ----
# y ~ N(mu_det, mu_det + mu_det^2/phi) on the observed cells; returns the same fields as cm_ekf_season
# so the two stages share every downstream consumer.
cm_det_season = function(y, f, S0, R0, c, b, phi, vax_day = NA, vax_frac = NULL, I0 = f$I0, engine = "R"){
  sim = if (engine == "R" || !exists("cm_simulate_season_engine", mode = "function")) cm_simulate_season(f, S0, R0, nrow(y), vax_day, vax_frac, I0)
        else cm_simulate_season_engine(f, S0, R0, nrow(y), vax_day, vax_frac, I0, engine = engine)
  mu = cm_mu(sim$inc, f, c, b); ok = is.finite(y)
  ll = sum(dnorm(y[ok], mu[ok], sqrt(mu[ok] + mu[ok]^2 / phi), log = TRUE))
  list(loglik = ll, mu_pred = mu, I = NULL, S = NULL)
}

# ---- |-penalised negative log-likelihood over a country's seasons ----
# stage = "det" scores the deterministic model (q unused), "ekf" the filter.
cm_negll = function(theta, cd, f, settings, R0_free = TRUE, return_fit = FALSE, stage = "ekf"){
  K = length(cd$seasons); p = cm_unpack(theta, K, R0_free, settings)
  # a wild line-search step (BFGS from a poor start) can overflow the exp/2^ transforms; return the
  # failure value rather than letting eigen() or the engines error out and kill the whole start
  if (!return_fit && !all(is.finite(c(p$S0, p$R0, p$I0, p$c, p$sigma, p$b, p$phi, p$q)))) return(1e10)
  engine = if (is.null(settings$engine) || !exists("cm_ekf_season_engine", mode = "function")) "R" else settings$engine
  f = cm_fixed_sigma(f, p$sigma)                                         # age-susceptibility profile (identity when off)
  ll = 0; filt = vector("list", K)
  for (s in seq_len(K)){
    bs = cm_b_season(p, cd, s)
    e = if (stage == "det") cm_det_season(cd$y[[s]], f, p$S0, p$R0[s], p$c[s, ], bs, p$phi, cd$vax_day, cd$vax_frac[[s]], I0 = p$I0[s], engine = engine)
        else if (engine == "R") cm_ekf_season(cd$y[[s]], f, p$S0, p$R0[s], p$c[s, ], bs, p$phi, p$q, cd$vax_day, cd$vax_frac[[s]], I0 = p$I0[s])
        else cm_ekf_season_engine(cd$y[[s]], f, p$S0, p$R0[s], p$c[s, ], bs, p$phi, p$q, cd$vax_day, cd$vax_frac[[s]], engine = engine, I0 = p$I0[s])
    ll = ll + e$loglik; filt[[s]] = e
  }
  lp = cm_logprior(theta, K, R0_free, settings)
  if (return_fit) return(list(negll = -(ll + lp), loglik = ll, logprior = lp, params = p, filt = filt))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
}

# ---- |-starting values (and multi-start jitter sds) in the layout of cm_unpack / cm_logprior ----
# The three functions must agree slot for slot: cm_theta_start builds the vector, cm_unpack reads it,
# cm_logprior penalises it. tests/testthat/test-comp-model-core.R checks the agreement for every
# switch combination (a missing slot silently shifts every later parameter -- the objective then
# reads phi off the end of the vector and returns the 1e10 failure value at every start).
cm_theta_start = function(cd, settings, R0_free = TRUE){
  K = length(cd$seasons)
  # data-driven starts: reporting from the observed peak (peak weekly incidence ~2% of a group), baseline from the floor
  rate_all = unlist(lapply(cd$rates, function(m) as.numeric(m)))
  peak_rate = max(rate_all, na.rm = TRUE); floor_rate = max(as.numeric(quantile(rate_all[rate_all > 0], 0.1, na.rm = TRUE)), 1e-3)
  c0 = min(max(peak_rate / settings$rate_per / 0.02, 1e-4), 1)     # mu = c N C -> c ~ (peak rate/1e5) / peak C
  # per-season seed start from the observed onset: at the starting growth rate r0 the default seed on
  # day 1 reaches onset (~1e-3 infected) around week 12; a season whose pooled rate first exceeds 10%
  # of its peak in week o_s is shifted by (o_s - 12) weeks -> log I0_s = log(I0) - r0 * 7 * (o_s - 12)
  by_season_I0 = isTRUE(settings$I0_by_season)
  r0 = settings$gamma_per_day * (settings$R0_reference * 0.8 - 1)
  onset = vapply(cd$rates, function(m){ tot = rowSums(sweep(m, 2, cd$N / sum(cd$N), "*"), na.rm = TRUE); tot[rowSums(is.finite(m)) == 0] = NA
    pk = max(tot, na.rm = TRUE); w = which(tot >= 0.1 * pk)[1]; if (is.na(w)) 12 else w }, numeric(1))
  logI0_start = log(settings$I0_fraction) - r0 * 7 * (onset - 12)
  by_age_c = isTRUE(settings$c_by_age); by_season_c = isTRUE(settings$c_by_season)
  by_age_susc = isTRUE(settings$susc_by_age); by_src_b = isTRUE(settings$b_by_source)
  base = c(qlogis(0.8), if (R0_free) rep(log(settings$R0_reference), K), if (by_season_I0) logI0_start,
           log(c0), if (by_age_c) c(0, 0), if (by_season_c) rep(0, K), if (by_age_susc) c(0, 0),
           if (by_src_b) c(log(max(floor_rate / 10, 1e-3)), log(floor_rate)) else log(floor_rate),
           log(15), log(0.1))
  names(base) = c("logit_S0", if (R0_free) paste0("log_R0_", cd$seasons), if (by_season_I0) paste0("log_I0_", cd$seasons),
                  "log_c", if (by_age_c) c("log2c_young", "log2c_elderly"), if (by_season_c) paste0("logc_dev_", cd$seasons),
                  if (by_age_susc) c("log2susc_young", "log2susc_elderly"),
                  if (by_src_b) c("log_b_RespiCompass", "log_b_ERVISS") else "log_b", "log_phi", "log_q")
  jit = c(0.6, if (R0_free) rep(0.03, K), if (by_season_I0) rep(0.7, K), 0.5, if (by_age_c) c(0.5, 0.5), if (by_season_c) rep(0.3, K),
          if (by_age_susc) c(0.5, 0.5), if (by_src_b) c(0.5, 0.5) else 0.5, 0.5, 0.5)
  list(base = base, jit = jit)
}

# ---- |-fit one country ----
# start: an optional NAMED vector (e.g. another fit's stage1$theta) whose slots override the
# data-driven start where the names match -- the WARM START for nested variants (a fit with an extra
# switch starts at the simpler fit's optimum with the new slots at 0, so it can only do better).
fit_comp_model = function(cd, settings, R0_free = TRUE, n_starts = settings$n_starts, seed = 1, verbose = TRUE,
                          cores = max(1, min(n_starts, parallel::detectCores() - 1)), start = NULL){
  f = cm_fixed(cd$Cn, cd$N, settings); K = length(cd$seasons)
  st = cm_theta_start(cd, settings, R0_free); base = st$base; jit = st$jit
  if (!is.null(start)){ common = intersect(names(start), names(base)); base[common] = start[common] }
  stage1 = NULL
  if (isTRUE(settings$two_stage)){
    # stage 1: the deterministic model, multi-start (cheap: no Jacobians, no filter)
    stage1 = .fit_multistart(base, jit, cm_negll,
                             list(cd = cd, f = f, settings = settings, R0_free = R0_free, stage = "det"),
                             n_starts, seed, sprintf("fit_comp_model(%s, deterministic)", cd$country),
                             cores = cores, maxit = settings$optim_maxit)
    if (verbose) cat(sprintf("   stage 1 (deterministic): negll %.1f conv=%d\n", stage1$value, stage1$convergence))
    base = stage1$par; if (is.null(settings$q_fixed)) base["log_q"] = settings$prior_logq["mean"]   # stage 2 starts here
    jit  = jit * 0.25                                                    # ... with only small jitter around it
  }
  best = .fit_multistart(base, jit, cm_negll,
                         list(cd = cd, f = f, settings = settings, R0_free = R0_free, stage = "ekf"),
                         if (isTRUE(settings$two_stage)) min(2, n_starts) else n_starts, seed,
                         sprintf("fit_comp_model(%s)", cd$country), cores = cores, maxit = settings$optim_maxit)
  fit = cm_negll(best$par, cd, f, settings, R0_free, return_fit = TRUE)
  fit_det = if (!is.null(stage1)) cm_negll(stage1$par, cd, f, settings, R0_free, return_fit = TRUE, stage = "det") else NULL
  # curvature-based standard errors on the transformed scale (Laplace)
  H = tryCatch(optimHess(best$par, cm_negll, cd = cd, f = f, settings = settings, R0_free = R0_free), error = function(e) NULL)
  se = if (!is.null(H)) tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(best$par))) else rep(NA_real_, length(best$par))
  names(se) = names(base)
  # deterministic curves + attack rates at the optimum (the mechanistic mean, no filtering)
  f_opt = cm_fixed_sigma(f, fit$params$sigma)
  det = lapply(seq_len(K), function(s){
    sim = cm_simulate_season(f_opt, fit$params$S0, fit$params$R0[s], cd$n_weeks[s], cd$vax_day, cd$vax_frac[[s]], I0 = fit$params$I0[s])
    list(mu = cm_mu(sim$inc, f_opt, fit$params$c[s, ], cm_b_season(fit$params, cd, s)), attack = sim$attack)
  })
  if (verbose) cat(sprintf("[%s] S0=%.3f  R0: %s  log10 I0: %s  c=%s b=%.2f phi=%.1f q=%.3f  negll=%.1f conv=%d\n",
                           cd$country, fit$params$S0, paste(sprintf("%.3f", fit$params$R0), collapse=" "),
                           paste(sprintf("%.1f", log10(fit$params$I0)), collapse=" "),
                           paste(sprintf("%.3g", fit$params$c_age), collapse="/"), mean(fit$params$b), fit$params$phi, fit$params$q, best$value, best$convergence))
  list(country = cd$country, seasons = cd$seasons, groups = cd$groups, N = cd$N, R0_free = R0_free,
       theta = best$par, se = se, params = fit$params, negll = best$value, loglik = fit$loglik,
       convergence = best$convergence, settings = settings,
       stage1 = if (!is.null(fit_det)) list(theta = stage1$par, params = fit_det$params, negll = stage1$value,
                                            mu = lapply(fit_det$filt, `[[`, "mu_pred")) else NULL,
       y = cd$y, rates = cd$rates, mu_filt = lapply(fit$filt, `[[`, "mu_pred"),
       I_filt = lapply(fit$filt, `[[`, "I"), S_filt = lapply(fit$filt, `[[`, "S"),
       mu_det = lapply(det, `[[`, "mu"), attack = unname(do.call(rbind, lapply(det, `[[`, "attack"))),
       R_eff = fit$params$R0 * fit$params$S0, vax = cd$vax, contact_source = cd$contact_source,
       source_by_season = cd$source_by_season)
}

# ---- |-tidy per-season summary of a fit (one row per season) ----
summarise_comp_fit = function(fit){
  K = length(fit$seasons)
  data.frame(country = fit$country, season = fit$seasons, S0 = fit$params$S0, R0 = fit$params$R0, I0 = fit$params$I0,
             R_eff = fit$R_eff, c = paste(signif(fit$params$c_age, 3), collapse = "/"), c_season = fit$params$c_season,
             susc = paste(signif(fit$params$sigma, 3), collapse = "/"),
             b = vapply(seq_len(K), function(s) cm_b_season(fit$params, fit, s), numeric(1)), phi = fit$params$phi, q = fit$params$q,
             attack_young = fit$attack[, 1], attack_medium = fit$attack[, 2], attack_elderly = fit$attack[, 3],
             cor = vapply(seq_len(K), function(s){ y = as.numeric(fit$y[[s]]); m = as.numeric(fit$mu_filt[[s]])
               ok = is.finite(y) & is.finite(m); if (sum(ok) > 2) cor(y[ok], m[ok]) else NA_real_ }, numeric(1)),
             vax_cov_65 = fit$vax$coverage, vax_provenance = fit$vax$provenance, stringsAsFactors = FALSE)
}
