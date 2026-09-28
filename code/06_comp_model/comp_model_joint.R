# comp_model_joint.R -- the JOINT fit: all countries and seasons at once, R0_s shared across countries
#
# ASSUMPTIONS.md C4/E2/I1 stage iv. Parameters:
#   log_R0[s]           one per SEASON, shared by every country (the Europe-wide transmissibility)
#   logit_S0[c]         one per COUNTRY, shared across its seasons
#   log_c[c], log_b[c], log_phi[c]   per-country reporting, baseline, dispersion
#   log_q               one process-noise scale for all
# Priors as in the per-country fit. The likelihood sums the EKF innovation log-likelihoods of every
# country-season; countries are independent GIVEN R0_s, so the per-country terms are evaluated in
# parallel (forked workers) with the C++ engine. Starting values come from the per-country fits
# (their S0, c, b, phi; R0_s from the across-country mean of the per-country R0_{c,s}).
#
# Cost (C++ engine, 4 cores): one joint likelihood for 5 eight-season countries ~0.1 s, a numerical
# gradient over ~30 parameters ~6 s, a fit ~30-60 min. For all 25 panel countries (~130 parameters)
# expect hours; analytical gradients would be the next step.

# layout: log_R0[K] | logit_S0[Cc] | log_c[Cc] | log_b[Cc] | log_phi[Cc] | log_q | [log_I0 per country-season, in cds order]
cm_joint_pack = function(fits, seasons_all, settings){
  K = length(seasons_all); Cc = length(fits)
  R0_by_season = sapply(seasons_all, function(s) mean(vapply(fits, function(f){ i = match(s, f$seasons); if (is.na(i)) NA_real_ else f$params$R0[i] }, numeric(1)), na.rm = TRUE))
  R0_by_season[!is.finite(R0_by_season)] = settings$R0_reference
  theta = c(log(R0_by_season),
            vapply(fits, function(f) qlogis(f$params$S0), numeric(1)),
            vapply(fits, function(f) log(f$params$c), numeric(1)),
            vapply(fits, function(f) log(f$params$b), numeric(1)),
            vapply(fits, function(f) log(f$params$phi), numeric(1)),
            log(mean(vapply(fits, function(f) f$params$q, numeric(1)))),
            if (isTRUE(settings$I0_by_season)) unlist(lapply(fits, function(f) log(f$params$I0))))
  names(theta) = c(paste0("log_R0_", seasons_all), paste0("logit_S0_", names(fits)), paste0("log_c_", names(fits)),
                   paste0("log_b_", names(fits)), paste0("log_phi_", names(fits)), "log_q",
                   if (isTRUE(settings$I0_by_season)) unlist(lapply(names(fits), function(cc) paste0("log_I0_", cc, "_", fits[[cc]]$seasons))))
  theta
}

cm_joint_unpack = function(theta, K, Cc, n_cs = NULL, settings = NULL){
  out = list(R0 = exp(theta[1:K]), S0 = plogis(theta[K + 1:Cc]), c = exp(theta[K + Cc + 1:Cc]),
             b = exp(theta[K + 2*Cc + 1:Cc]), phi = exp(theta[K + 3*Cc + 1:Cc]), q = exp(theta[K + 4*Cc + 1]))
  if (!is.null(settings) && isTRUE(settings$I0_by_season)){               # per country-season seeds, split by country
    v = exp(theta[K + 4*Cc + 1 + seq_len(sum(n_cs))]); out$I0 = split(v, rep(seq_len(Cc), n_cs))
  } else out$I0 = lapply(n_cs, function(n) rep(settings$I0_fraction, n))
  out
}

cm_joint_logprior = function(theta, K, Cc, settings, n_cs = NULL){
  lp = sum(dnorm(theta[1:K], log(settings$R0_reference), settings$prior_logR0_sd, log = TRUE)) +
  sum(dnorm(theta[K + 1:Cc], settings$prior_logitS0["mean"], settings$prior_logitS0["sd"], log = TRUE)) +
  sum(dnorm(theta[K + 3*Cc + 1:Cc], settings$prior_logphi["mean"], settings$prior_logphi["sd"], log = TRUE)) +
  dnorm(theta[K + 4*Cc + 1], settings$prior_logq["mean"], settings$prior_logq["sd"], log = TRUE)
  if (isTRUE(settings$I0_by_season)) lp = lp + sum(dnorm(theta[K + 4*Cc + 1 + seq_len(sum(n_cs))], log(settings$I0_fraction), settings$prior_logI0_sd, log = TRUE))
  lp
}

# cds: named list of build_comp_data() outputs; seasons_all: the union of seasons (defines R0_s slots)
cm_joint_negll = function(theta, cds, fixed, seasons_all, settings, engine = "cpp", cores = 1, return_fit = FALSE){
  K = length(seasons_all); Cc = length(cds); n_cs = vapply(cds, function(cd) length(cd$seasons), integer(1))
  p = cm_joint_unpack(theta, K, Cc, n_cs, settings)
  one_country = function(i){
    cd = cds[[i]]; f = fixed[[i]]; ll = 0; filt = vector("list", length(cd$seasons))
    for (s in seq_along(cd$seasons)){
      k = match(cd$seasons[s], seasons_all)
      e = cm_ekf_season_engine(cd$y[[s]], f, p$S0[i], p$R0[k], p$c[i], p$b[i], p$phi[i], p$q, cd$vax_day, cd$vax_frac[[s]], engine = engine, I0 = p$I0[[i]][s])
      ll = ll + e$loglik; if (return_fit) filt[[s]] = e
    }
    list(ll = ll, filt = filt)
  }
  res = if (cores > 1) parallel::mclapply(seq_len(Cc), one_country, mc.cores = cores) else lapply(seq_len(Cc), one_country)
  ll = sum(vapply(res, function(r) if (inherits(r, "try-error") || is.null(r)) -1e10 else r$ll, numeric(1)))
  lp = cm_joint_logprior(theta, K, Cc, settings, n_cs)
  if (return_fit) return(list(negll = -(ll + lp), loglik = ll, params = p, filt = lapply(res, `[[`, "filt")))
  if (!is.finite(ll + lp)) return(1e10)
  -(ll + lp)
}

fit_comp_model_joint = function(cds, fits, settings, engine = "cpp", cores = max(1, parallel::detectCores() - 1),
                                maxit = 400, verbose = TRUE){
  if (engine == "cpp") cm_load_cpp()
  seasons_all = sort(unique(unlist(lapply(cds, `[[`, "seasons"))))
  fixed = lapply(cds, function(cd) cm_fixed(cd$Cn, cd$N, settings))
  theta0 = cm_joint_pack(fits[names(cds)], seasons_all, settings)
  t0 = Sys.time()
  opt = optim(theta0, cm_joint_negll, cds = cds, fixed = fixed, seasons_all = seasons_all, settings = settings,
              engine = engine, cores = cores, method = "BFGS", control = list(maxit = maxit, reltol = 1e-8))
  fit = cm_joint_negll(opt$par, cds, fixed, seasons_all, settings, engine, cores, return_fit = TRUE)
  H = tryCatch(optimHess(opt$par, cm_joint_negll, cds = cds, fixed = fixed, seasons_all = seasons_all, settings = settings,
                         engine = engine, cores = cores), error = function(e) NULL)
  se = if (!is.null(H)) tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(opt$par))) else rep(NA_real_, length(opt$par))
  names(se) = names(theta0)
  if (verbose) cat(sprintf("joint fit: %d countries, %d seasons, %d parameters, negll %.1f -> %.1f, conv=%d, %.0f min\n",
                           length(cds), length(seasons_all), length(theta0), cm_joint_negll(theta0, cds, fixed, seasons_all, settings, engine, cores),
                           opt$value, opt$convergence, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(countries = names(cds), seasons = seasons_all, theta = opt$par, se = se, params = fit$params,
       negll = opt$value, loglik = fit$loglik, convergence = opt$convergence, filt = fit$filt, settings = settings)
}

# ---- |-tidy tables from a joint fit ----
summarise_joint_fit = function(jf){
  K = length(jf$seasons); Cc = length(jf$countries)
  list(seasons = data.frame(season = jf$seasons, R0 = jf$params$R0,
                            lo = exp(jf$theta[1:K] - 1.96 * jf$se[1:K]), hi = exp(jf$theta[1:K] + 1.96 * jf$se[1:K])),
       countries = data.frame(country = jf$countries, S0 = jf$params$S0,
                              lo = plogis(jf$theta[K + 1:Cc] - 1.96 * jf$se[K + 1:Cc]), hi = plogis(jf$theta[K + 1:Cc] + 1.96 * jf$se[K + 1:Cc]),
                              c = jf$params$c, b = jf$params$b, phi = jf$params$phi))
}
