# joint_model.R -- data assembly, parameter packing, and the fitting strategy for the joint model.
# The model itself is MODEL.md; the log-posterior is joint_model.cpp; this file is the driver.
#
# Fitting strategy, and why. The likelihood is SEPARABLE: only R0_s, the season deviations and the one
# elderly susceptibility appear in more than one country's terms. Everything else is local to a single
# country. A flat finite-difference gradient over ~184 parameters costs 185 full-likelihood
# evaluations; a block sweep costs about 33, because a country's own likelihood is a twelfth of the
# joint one. So we alternate: every country's local block optimised in parallel with the shared block
# held fixed, then the shared block with the locals held fixed, to convergence, then one joint polish.

suppressMessages({library(Rcpp); library(dplyr)})

# ---- |-settings: the fixed constants and the priors, in one place ----
jm_settings = function(){
  list(
    # fixed, not fitted (MODEL.md)
    gamma_per_day = 0.2777778,     # 3.6-day mean infectious period, from the project notes
    ve_inf = 0.25, ve_ili = 0.20, ve_spread = 0.20,
    rate_per = 1e5, vax_day = 62L, # the 65+ vaccination pulse, 1 October
    # priors, all on the unconstrained scale the fit works in
    pr_R0_mean = log(1.5), pr_R0_sd = 0.15,   # widened from the pilot's 0.05: R0_s is now shared, so
                                              # 86 waves inform 8 numbers and the data can carry it
    pr_S0_mean = qlogis(0.75), pr_S0_sd = 1,
    pr_sigma_mean = 1, pr_sigma_sd = 0.5,     # log2 sigma_eld ~ N(1, .5): centre 2x, 95% band 1.0-4.0x
    pr_delta_sd = 0.5,                        # season observation deviation, constrained to average 0
    pr_off_sd = 1,                            # log2 age reporting offsets
    pr_c_mean = log(0.05), pr_c_sd = 2,       # weak; closes the 'no epidemic' trap at c -> 0
    # The data's own week-to-week scatter around a 3-week mean implies phi ~ 3.7 (median over 257
    # country-season-age series), so that is the centre. The sd is deliberately WIDE: without the
    # filter the mean is rigid, so phi must be free to absorb model misfit and REPORT it. A tight
    # prior would force phi up and push the misfit into R0 and S0 instead, which is worse. The gap
    # between the fitted phi and 3.7 is therefore the model-adequacy diagnostic (jm_adequacy).
    pr_phi_mean = log(4), pr_phi_sd = 1.0,
    pr_b_mean = log(1), pr_b_sd = 2,          # NEW: the pilot left b unpenalised, which left a
                                              # single-season baseline effectively free
    pr_I0_mean = log(10^-6.5), pr_I0_sd = 3
  )
}

jm_load_cpp = function(dir = "output/joint_model/cpp_cache"){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(here::here("code/07_joint_model/joint_model.cpp"), cacheDir = dir, verbose = FALSE)
  invisible(TRUE)
}

# ---- |-assemble every country into the flat structure the C++ expects ----
# Reuses the compartmental pilot's data layer (build_comp_data): the panel stitch, the 3-group
# collapse of the contact matrix, the populations and the 65+ coverage are already settled there.
# Counts are ROUNDED here: the panel holds rates, the count scale is a device (MODEL.md), and the
# negative binomial wants integers.
jm_build_data = function(countries, models_in, demo, set = jm_settings(), min_seasons = 6L, verbose = TRUE){
  cds = list()
  for (cc in countries){
    cd = tryCatch(build_comp_data(cc, models_in, demo, comp_model_settings()), error = function(e) NULL)
    if (is.null(cd)){ if (verbose) cat("  skip", cc, "(no data)\n"); next }
    if (length(cd$seasons) < min_seasons){ if (verbose) cat("  skip", cc, sprintf("(%d seasons < %d)\n", length(cd$seasons), min_seasons)); next }
    cds[[cc]] = cd
  }
  stopifnot(length(cds) > 0)
  seasons_all = sort(unique(unlist(lapply(cds, `[[`, "seasons"))))
  C = length(cds); S = length(seasons_all)

  # per country: the data sources it actually has. The pilot created a slot for both sources
  # unconditionally, which left an inert parameter wherever a source was absent (ES, NO) and made the
  # Hessian exactly singular. Here the slots follow the data.
  srcs = lapply(cds, function(cd) sort(unique(unname(cd$source_by_season))))
  n_src = vapply(srcs, length, integer(1))
  n_cs_of_country = vapply(cds, function(cd) length(cd$seasons), integer(1))
  n_local = 5L + n_src + n_cs_of_country                 # S0, c, 2 offsets, phi, baselines, seeds
  n_shared = 2L * S                                      # S log-R0, S-1 free deviations, 1 log2 sigma
  off_country = as.integer(cumsum(c(n_shared, head(n_local, -1))))   # 0-based starts
  n_par = as.integer(n_shared + sum(n_local))

  y = list(); cs_country = integer(0); cs_season = integer(0); cs_src = integer(0)
  cs_pos = integer(0); vax_eld = numeric(0); lgamma_y1 = numeric(0)
  cs_label = character(0); n_weeks = integer(0)
  for (ic in seq_along(cds)){
    cd = cds[[ic]]
    for (k in seq_along(cd$seasons)){
      ym = round(cd$y[[k]]); ym[ym < 0] = 0
      y[[length(y) + 1]] = ym
      cs_country = c(cs_country, ic - 1L)
      cs_season  = c(cs_season, match(cd$seasons[k], seasons_all) - 1L)
      cs_src     = c(cs_src, match(unname(cd$source_by_season[k]), srcs[[ic]]) - 1L)
      cs_pos     = c(cs_pos, k - 1L)
      vax_eld    = c(vax_eld, cd$vax$coverage[k])
      lgamma_y1  = c(lgamma_y1, sum(lgamma(ym[is.finite(ym)] + 1)))
      cs_label   = c(cs_label, paste(names(cds)[ic], cd$seasons[k]))
      n_weeks    = c(n_weeks, nrow(ym))
    }
  }
  cs_of_country = lapply(seq_along(cds), function(ic) as.integer(which(cs_country == ic - 1L) - 1L))

  d = c(list(
    n_country = C, n_season = S, n_cs = length(y), n_par = n_par,
    countries = names(cds), seasons = seasons_all, sources = srcs,
    Cn = lapply(cds, `[[`, "Cn"), N = lapply(cds, `[[`, "N"), groups = cds[[1]]$groups,
    n_src = as.integer(n_src), n_cs_of_country = as.integer(n_cs_of_country),
    n_local = as.integer(n_local), off_country = off_country,
    cs_of_country = cs_of_country, cs_country = cs_country, cs_season = cs_season,
    cs_src = cs_src, cs_pos = cs_pos, cs_label = cs_label, n_weeks = as.integer(n_weeks),
    vax_eld = vax_eld, lgamma_y1 = lgamma_y1, y = y,
    rates = unlist(lapply(cds, `[[`, "rates"), recursive = FALSE),
    gamma = set$gamma_per_day, ve_inf = set$ve_inf, ve_ili = set$ve_ili,
    ve_spread = set$ve_spread, rate_per = set$rate_per, vax_day = set$vax_day
  ), set[grep("^pr_", names(set))])
  if (verbose) cat(sprintf("%d countries, %d seasons, %d country-seasons, %d parameters (%d shared, %d local)\n",
                           C, S, length(y), n_par, n_shared, sum(n_local)))
  d
}

# ---- |-the parameter vector: names, index blocks, starting values ----
jm_par_names = function(d){
  S = d$n_season
  nm = c(paste0("log_R0_", d$seasons), paste0("delta_", head(d$seasons, S - 1)), "log2_sigma_eld")
  for (ic in seq_len(d$n_country)){
    cc = d$countries[ic]
    nm = c(nm, paste0(cc, c(":logit_S0", ":log_c", ":off_young", ":off_eld", ":log_phi")),
           paste0(cc, ":log_b_", d$sources[[ic]]),
           paste0(cc, ":log_I0_", d$seasons[d$cs_season[d$cs_of_country[[ic]] + 1L] + 1L]))
  }
  nm
}
jm_blocks = function(d){
  S = d$n_season
  list(shared = seq_len(2L * S),
       local = lapply(seq_len(d$n_country), function(ic) d$off_country[ic] + seq_len(d$n_local[ic])))
}

jm_theta0 = function(d, set = jm_settings()){
  th = numeric(d$n_par)
  S = d$n_season
  th[seq_len(S)] = set$pr_R0_mean                       # all seasons at 1.5
  th[S + seq_len(S - 1)] = 0                            # no season deviation
  th[2L * S] = set$pr_sigma_mean                        # elderly susceptibility at 2x
  r0 = set$gamma_per_day * (1.5 * 0.8 - 1)              # a plausible early growth rate, per day
  for (ic in seq_len(d$n_country)){
    base = d$off_country[ic]
    ics = d$cs_of_country[[ic]] + 1L
    rate_all = unlist(lapply(ics, function(i) as.numeric(d$rates[[i]])))
    peak = max(rate_all, na.rm = TRUE)
    floor_r = max(as.numeric(quantile(rate_all[rate_all > 0], 0.10, na.rm = TRUE)), 1e-3)
    c0 = min(max(peak / d$rate_per / 0.02, 1e-4), 1)
    th[base + 1] = qlogis(0.80)
    th[base + 2] = log(c0)
    th[base + 3] = 0; th[base + 4] = 0
    th[base + 5] = set$pr_phi_mean
    th[base + 5 + seq_len(d$n_src[ic])] = log(floor_r)
    # the seed start comes from the OBSERVED onset: at growth rate r0 the default seed reaches onset
    # around week 12, so a season whose pooled rate first passes 10% of its peak in week o is shifted
    # by (o - 12) weeks. This is the only handle on arrival time, so a good start matters.
    onset = vapply(ics, function(i){
      m = d$rates[[i]]; tot = rowSums(sweep(m, 2, d$N[[ic]] / sum(d$N[[ic]]), "*"), na.rm = TRUE)
      tot[rowSums(is.finite(m)) == 0] = NA
      w = which(tot >= 0.1 * max(tot, na.rm = TRUE))[1]; if (is.na(w)) 12 else w }, numeric(1))
    th[base + 5 + d$n_src[ic] + seq_along(ics)] = set$pr_I0_mean - r0 * 7 * (onset - 12)
  }
  names(th) = jm_par_names(d)
  th
}

# ---- |-fit by block coordinate descent, then one joint polish ----
jm_fit = function(d, theta0 = jm_theta0(d), max_sweeps = 15L, tol = 0.05, maxit_local = 400L,
                  maxit_shared = 400L, polish_maxit = 300L, cores = max(1L, parallel::detectCores() - 1L),
                  verbose = TRUE){
  bl = jm_blocks(d); th = theta0
  t_start = Sys.time()
  obj = jm_negll_cpp(th, d)
  trace = data.frame(sweep = 0L, step = "start", negll = obj, seconds = 0)
  if (verbose) cat(sprintf("start: negll %.2f\n", obj))
  conv_local = integer(0)
  for (sw in seq_len(max_sweeps)){
    prev = obj
    # --- every country's local block, with the shared block fixed. Independent, so run in parallel.
    fits = parallel::mclapply(seq_len(d$n_country), function(ic){
      idx = bl$local[[ic]]
      f = function(x){ t2 = th; t2[idx] = x; jm_country_negll_cpp(t2, d, ic - 1L) }
      o = tryCatch(optim(th[idx], f, method = "BFGS", control = list(maxit = maxit_local, reltol = 1e-10)),
                   error = function(e) NULL)
      if (is.null(o)) NULL else list(par = o$par, value = o$value, conv = o$convergence)
    }, mc.cores = min(cores, d$n_country))
    for (ic in seq_len(d$n_country)) if (!is.null(fits[[ic]])) th[bl$local[[ic]]] = fits[[ic]]$par
    conv_local = vapply(fits, function(f) if (is.null(f)) 99L else as.integer(f$conv), integer(1))
    obj_l = jm_negll_cpp(th, d)
    trace = rbind(trace, data.frame(sweep = sw, step = "local", negll = obj_l,
                                    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    # --- the shared block, with every local block fixed
    idx = bl$shared
    fs = function(x){ t2 = th; t2[idx] = x; jm_negll_cpp(t2, d) }
    os = tryCatch(optim(th[idx], fs, method = "BFGS", control = list(maxit = maxit_shared, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(os)) th[idx] = os$par
    obj = jm_negll_cpp(th, d)
    trace = rbind(trace, data.frame(sweep = sw, step = "shared", negll = obj,
                                    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    if (verbose) cat(sprintf("sweep %2d: after locals %.2f, after shared %.2f  (gain %.2f)\n",
                             sw, obj_l, obj, prev - obj))
    if (prev - obj < tol) break
  }
  # --- one joint polish: all parameters together, to clear any residual cross-block curvature
  op = tryCatch(optim(th, function(x) jm_negll_cpp(x, d), method = "BFGS",
                      control = list(maxit = polish_maxit, reltol = 1e-11)), error = function(e) NULL)
  polish_gain = NA_real_
  if (!is.null(op) && is.finite(op$value) && op$value < obj){ polish_gain = obj - op$value; th = op$par; obj = op$value }
  names(th) = jm_par_names(d)
  trace = rbind(trace, data.frame(sweep = 99L, step = "polish", negll = obj,
                                  seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  secs = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) cat(sprintf("done in %.1f s: negll %.2f (polish gained %.2f), loglik %.2f\n",
                           secs, obj, polish_gain, jm_loglik_cpp(th, d)))
  list(theta = th, negll = obj, loglik = jm_loglik_cpp(th, d), seconds = secs, trace = trace,
       sweeps = if (exists("sw")) sw else NA_integer_, conv_local = conv_local,
       conv_shared = if (is.null(os)) 99L else os$convergence,
       conv_polish = if (is.null(op)) 99L else op$convergence, polish_gain = polish_gain, d = d)
}

# ---- |-tidy the parameters ----
jm_unpack = function(th, d){
  S = d$n_season
  dev_free = th[S + seq_len(S - 1)]
  out = list(R0 = setNames(exp(th[seq_len(S)]), d$seasons),
             delta = setNames(c(dev_free, -sum(dev_free)), d$seasons),
             sigma_eld = unname(2^th[2L * S]))
  cty = lapply(seq_len(d$n_country), function(ic){
    base = d$off_country[ic]; ics = d$cs_of_country[[ic]] + 1L
    list(S0 = plogis(th[base + 1]), c = exp(th[base + 2]),
         off_young = th[base + 3], off_eld = th[base + 4], phi = exp(th[base + 5]),
         b = setNames(exp(th[base + 5 + seq_len(d$n_src[ic])]), d$sources[[ic]]),
         I0 = setNames(exp(th[base + 5 + d$n_src[ic] + seq_along(ics)]), d$seasons[d$cs_season[ics] + 1L]))
  })
  names(cty) = d$countries
  out$country = cty
  out
}

jm_summary_country = function(fit){
  p = jm_unpack(fit$theta, fit$d); d = fit$d
  data.frame(country = d$countries,
             S0 = vapply(p$country, `[[`, numeric(1), "S0"),
             c_adult = vapply(p$country, `[[`, numeric(1), "c"),
             rel_young = 2^vapply(p$country, `[[`, numeric(1), "off_young"),
             rel_elderly = 2^vapply(p$country, `[[`, numeric(1), "off_eld"),
             phi = vapply(p$country, `[[`, numeric(1), "phi"),
             n_seasons = d$n_cs_of_country, row.names = NULL)
}
jm_summary_season = function(fit){
  p = jm_unpack(fit$theta, fit$d)
  data.frame(season = fit$d$seasons, R0 = unname(p$R0), deviation = unname(p$delta),
             reporting_mult = exp(unname(p$delta)), row.names = NULL)
}

# ---- |-base-R reference implementation of the SAME model, for the identity test ----
# Deliberately written straight from the maths rather than by translating the C++, and using R's own
# eigen() for the spectral radius rather than the C++ power iteration, so the two agree only if both
# are right. tests/testthat/test-joint-model.R requires 1e-10.
jm_negll_R = function(th, d){
  S = d$n_season; A = 3L
  R0 = exp(th[seq_len(S)])
  dev_free = th[S + seq_len(S - 1)]
  dev = exp(c(dev_free, -sum(dev_free)))
  sigma = c(1, 1, 2^th[2L * S])
  dn = function(x, m, s) -0.5 * ((x - m) / s)^2 - log(s) - 0.5 * log(2 * pi)
  lp = sum(dn(th[seq_len(S)], d$pr_R0_mean, d$pr_R0_sd)) +
       sum(dn(c(dev_free, -sum(dev_free)), 0, d$pr_delta_sd)) +
       dn(th[2L * S], d$pr_sigma_mean, d$pr_sigma_sd)
  for (ic in seq_len(d$n_country)){
    base = d$off_country[ic]; nsrc = d$n_src[ic]
    S0 = plogis(th[base + 1]); cc = exp(th[base + 2])
    c_age = cc * 2^c(th[base + 3], 0, th[base + 4])
    phi = exp(th[base + 5])
    b = exp(th[base + 5 + seq_len(nsrc)])
    I0v = exp(th[base + 5 + nsrc + seq_len(d$n_cs_of_country[ic])])
    Cs = sweep(d$Cn[[ic]], 1, sigma, "*")
    Cs = Cs / max(abs(eigen(Cs, only.values = TRUE)$values))
    N = d$N[[ic]]
    for (ics in d$cs_of_country[[ic]] + 1L){
      y = d$y[[ics]]; nw = nrow(y); s = d$cs_season[ics] + 1L
      beta = R0[s] * d$gamma
      Su = rep(S0, A); Iu = rep(I0v[d$cs_pos[ics] + 1L], A); Sv = rep(0, A); Iv = rep(0, A)
      vax = c(0, 0, d$vax_eld[ics]); inc = matrix(0, nw, A); day = 0L
      for (t in seq_len(nw)){
        acc = rep(0, A)
        for (k in 1:7){
          day = day + 1L
          lam = beta * as.numeric(Cs %*% (Iu + (1 - d$ve_spread) * Iv))
          nu = lam * Su; nv = (1 - d$ve_inf) * lam * Sv
          Su = Su - nu; Iu = Iu + nu - d$gamma * Iu
          Sv = Sv - nv; Iv = Iv + nv - d$gamma * Iv
          acc = acc + nu + (1 - d$ve_ili) * nv
          if (day == d$vax_day){ mv = vax * Su; Su = Su - mv; Sv = Sv + mv }
          Su = pmin(pmax(Su, 0), 1); Iu = pmin(pmax(Iu, 0), 1)
          Sv = pmin(pmax(Sv, 0), 1); Iv = pmin(pmax(Iv, 0), 1)
        }
        inc[t, ] = acc
      }
      mu = sweep(inc, 2, c_age * dev[s] * N, "*") +
           matrix(b[d$cs_src[ics] + 1L] * N / d$rate_per, nw, A, byrow = TRUE)
      mu[mu < 1e-10] = 1e-10
      ok = is.finite(y)
      lp = lp + sum(lgamma(y[ok] + phi) - lgamma(phi) + phi * (log(phi) - log(phi + mu[ok])) +
                    y[ok] * (log(mu[ok]) - log(phi + mu[ok]))) - d$lgamma_y1[ics]
      lp = lp + dn(log(I0v[d$cs_pos[ics] + 1L]), d$pr_I0_mean, d$pr_I0_sd)
    }
    lp = lp + dn(th[base + 1], d$pr_S0_mean, d$pr_S0_sd) + dn(th[base + 2], d$pr_c_mean, d$pr_c_sd) +
         dn(th[base + 3], 0, d$pr_off_sd) + dn(th[base + 4], 0, d$pr_off_sd) +
         dn(th[base + 5], d$pr_phi_mean, d$pr_phi_sd) +
         sum(dn(th[base + 5 + seq_len(nsrc)], d$pr_b_mean, d$pr_b_sd))
  }
  -lp
}

# ---- |-convergence and identifiability indicators ----
# The likelihood-only Hessian is the curvature of the fit-quality landscape with the priors switched
# off: steep directions are combinations the data determine, near-flat ones are combinations the data
# cannot see. Comparing it with the penalised Hessian shows what the priors are holding up.
# `contraction` = 1 - sd_posterior/sd_prior: 0 means the posterior is just the prior, 1 means the data
# did the work. `sd_data` is the data-only precision with the other parameters held fixed.
jm_identifiability = function(fit, verbose = TRUE){
  d = fit$d; th = fit$theta; nm = names(th)
  if (verbose) cat("Hessians (", length(th), "parameters, ~", round(2 * length(th)^2 / 1000), "k evaluations )...\n")
  t0 = Sys.time()
  H_post = optimHess(th, function(x) jm_negll_cpp(x, d))
  H_lik  = optimHess(th, function(x) -jm_loglik_cpp(x, d))
  H_post = (H_post + t(H_post)) / 2; H_lik = (H_lik + t(H_lik)) / 2
  secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  S = d$n_season
  fam = ifelse(grepl("log_R0", nm), "R0 (season, shared)",
        ifelse(grepl("^delta_", nm), "season deviation (shared)",
        ifelse(grepl("log2_sigma", nm), "elderly susceptibility (global)",
        ifelse(grepl(":logit_S0", nm), "S0 (country)",
        ifelse(grepl(":log_c$", nm), "reporting c (country)",
        ifelse(grepl(":off_", nm), "age reporting offset",
        ifelse(grepl(":log_phi", nm), "dispersion phi",
        ifelse(grepl(":log_b_", nm), "baseline b",
        ifelse(grepl(":log_I0_", nm), "seed I0 (country-season)", nm)))))))))
  pr_sd = ifelse(grepl("log_R0", nm), d$pr_R0_sd,
          ifelse(grepl("^delta_", nm), d$pr_delta_sd,
          ifelse(grepl("log2_sigma", nm), d$pr_sigma_sd,
          ifelse(grepl(":logit_S0", nm), d$pr_S0_sd,
          ifelse(grepl(":log_c$", nm), d$pr_c_sd,
          ifelse(grepl(":off_", nm), d$pr_off_sd,
          ifelse(grepl(":log_phi", nm), d$pr_phi_sd,
          ifelse(grepl(":log_b_", nm), d$pr_b_sd,
          ifelse(grepl(":log_I0_", nm), d$pr_I0_sd, NA_real_)))))))))
  e_lik = eigen(H_lik, symmetric = TRUE, only.values = TRUE)$values
  e_post = eigen(H_post, symmetric = TRUE, only.values = TRUE)$values
  V_post = tryCatch(solve(H_post), error = function(e) NULL)
  sd_post = if (is.null(V_post)) rep(NA_real_, length(th)) else sqrt(pmax(diag(V_post), 0))
  sd_data = 1 / sqrt(pmax(diag(H_lik), 1e-300))
  # the shared block on its own, conditional on every local parameter: what pooling actually delivers
  sh = seq_len(2L * S)
  V_sh = tryCatch(solve(H_lik[sh, sh, drop = FALSE]), error = function(e) NULL)
  sd_shared_cond = if (is.null(V_sh)) rep(NA_real_, length(sh)) else sqrt(pmax(diag(V_sh), 0))
  tab = data.frame(parameter = nm, family = fam, prior_sd = pr_sd, sd_data = sd_data,
                   sd_post = sd_post, contraction = 1 - sd_post / pr_sd, row.names = NULL)
  famtab = tab %>% group_by(family) %>%
    summarise(n = n(), prior_sd = first(prior_sd), sd_data_med = median(sd_data),
              sd_post_med = median(sd_post), contraction_med = median(contraction),
              contraction_min = min(contraction), .groups = "drop") %>% arrange(desc(contraction_med))
  list(seconds = secs, eig_lik = e_lik, eig_post = e_post,
       flat_ratio = min(abs(e_lik)) / max(abs(e_lik)),
       n_near_flat = sum(abs(e_lik) < 1e-8 * max(abs(e_lik))),
       n_negative = sum(e_lik < -1e-6 * max(abs(e_lik))),
       post_pd = all(e_post > 0), table = tab, family = as.data.frame(famtab),
       sd_shared_cond = setNames(sd_shared_cond, nm[sh]))
}

# ---- |-model adequacy: is the fitted dispersion explainable as MEASUREMENT noise? ----
# The negative binomial has Var = mu + mu^2/phi, so at large counts the coefficient of variation is
# 1/sqrt(phi). The data's own scatter around a 3-week moving average is a lower bound on what any
# mean could achieve. If the fitted phi implies much more noise than that, the excess is the
# deterministic mean misfitting, not measurement -- and that is the case for adding the filter back.
jm_adequacy = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d)
  ma3 = function(v){ n = length(v); o = rep(NA_real_, n)
    if (n >= 3) for (i in 2:(n - 1)){ w = v[(i - 1):(i + 1)]; if (all(is.finite(w))) o[i] = mean(w) }; o }
  rows = lapply(seq_len(d$n_country), function(ic){
    cvs = unlist(lapply(d$cs_of_country[[ic]] + 1L, function(i){
      y = d$y[[i]]; vapply(seq_len(ncol(y)), function(a){
        v = y[, a]; m = ma3(v); ok = is.finite(v) & is.finite(m) & m > 20
        if (sum(ok) < 5) NA_real_ else sd(v[ok] / m[ok]) / sqrt(2/3) }, numeric(1)) }))
    cv_data = median(cvs, na.rm = TRUE)
    data.frame(country = d$countries[ic], phi_fitted = p$country[[ic]]$phi,
               cv_fitted = 1 / sqrt(p$country[[ic]]$phi), cv_data = cv_data,
               phi_data = 1 / cv_data^2, excess = (1 / sqrt(p$country[[ic]]$phi)) / cv_data)
  })
  do.call(rbind, rows)
}

# ---- |-observed vs fitted, as rates per 100 000, tidy ----
jm_tidy_fit = function(fit){
  d = fit$d; f = jm_fitted_cpp(fit$theta, d)
  do.call(rbind, lapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; N = d$N[[ic]]; nw = d$n_weeks[i]
    per100k = d$rate_per / N
    obs = sweep(d$y[[i]], 2, per100k, "*"); fitv = sweep(f$mu[[i]], 2, per100k, "*")
    data.frame(country = d$countries[ic], season = d$seasons[d$cs_season[i] + 1L],
               week = rep(seq_len(nw), 2 * length(N)),
               group = rep(rep(d$groups, each = nw), 2),
               what = rep(c("observed", "model"), each = nw * length(N)),
               value = c(as.numeric(obs), as.numeric(fitv)), stringsAsFactors = FALSE)
  }))
}
