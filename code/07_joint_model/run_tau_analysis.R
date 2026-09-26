# run_tau_analysis.R -- the spatial spread tau: is it there, is it by country, and what can the data
# tell apart among the three things that make a wave fat (R0, S0 and tau)? MODEL.md, "Spatial spread".
#
# Every step caches its result in output/joint_model/tau/ and is skipped when the file exists, so the
# analysis can be re-run piecewise (delete a file to redo that step). Full run: about two hours on four
# cores, most of it the by-country fits (6 min each) and the Hessians.
#
#   1. fits    tau pinned at 0 (a needle prior; must reproduce the tau-free optimum), shared tau (two
#              starts), tau by country (five starts: the landscape has more than one optimum), and
#              shared tau with R0 pinned at 1.7 (is tau held down by the S0 <= 1 ceiling?)
#   2. profile the shared tau held on a grid, everything else refitted
#   3. profile each country's tau held on a grid, that country's block refitted (conditional on the
#              shared block), recording its S0 as it goes: the tau-S0 trade-off, country by country
#   4. Hessians identifiability of the shared and the by-country fits
#   5. recovery a known shared tau (fitted, 7 and 14 days); a by-country world where big countries
#              spread (14 d) and small ones do not (3 d), fitted both ways
#   6. compare wave-level inference between the variants; figure 18
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_recovery.R")
source("code/07_joint_model/joint_compare.R"); jm_load_cpp()
suppressMessages({library(ggplot2); library(patchwork)})
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj
options(width = 200)
out = "output/joint_model/tau"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
cands = c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")
cores = max(1L, parallel::detectCores() - 1L)
cached = function(file, expr){ f = file.path(out, file); if (file.exists(f)) return(readRDS(f)); v = expr; saveRDS(v, f); v }
best_of = function(...){ fs = list(...); fs[[which.min(vapply(fs, `[[`, 0, "negll"))]] }
set_c = modifyList(jm_settings(), list(tau_by_country = TRUE))
taus_of = function(fit) vapply(jm_unpack(fit$theta, fit$d)$country, `[[`, 0, "tau")
set_country_tau = function(th, d, v){ for (ic in seq_len(d$n_country)) th[d$off_country[ic] + d$n_local[ic]] = v; th }
# shared layout -> by-country layout: drop the shared slot, give every country the same tau
to_country = function(th_s, d_s, d_c, log_tau){
  S = d_s$n_season; th = numeric(d_c$n_par); th[seq_len(2 * S - 1)] = th_s[seq_len(2 * S - 1)]
  for (ic in seq_len(d_c$n_country))
    th[d_c$off_country[ic] + seq_len(d_c$n_local[ic])] = c(th_s[d_s$off_country[ic] + seq_len(d_s$n_local[ic])], log_tau)
  names(th) = jm_par_names(d_c); th
}

# ---- 1. fits ----
f0 = cached("fit_tau0.rds", {
  s = modifyList(jm_settings(), list(pr_tau_mean = -30, pr_tau_sd = 0.01))
  d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  jm_fit(d, theta0 = jm_theta0(d, s), cores = cores, verbose = FALSE) })
dS = jm_build_data(cands, models_in, demo, verbose = FALSE)
fS = cached("fit_tau_shared.rds", {
  warm = f0$theta; warm[["log_tau"]] = log(7)
  best_of(jm_fit(dS, cores = cores, verbose = FALSE), jm_fit(dS, theta0 = warm, cores = cores, verbose = FALSE)) })
dC = jm_build_data(cands, models_in, demo, set = set_c, verbose = FALSE)
fC = cached("fit_tau_country.rds", {
  a = jm_fit(dC, theta0 = jm_theta0(dC, set_c), cores = cores, verbose = FALSE)
  b = jm_fit(dC, theta0 = to_country(fS$theta, dS, dC, fS$theta[["log_tau"]]), cores = cores, verbose = FALSE)
  best = best_of(a, b)
  c3 = jm_fit(dC, theta0 = best$theta, cores = cores, verbose = FALSE)
  c4 = jm_fit(dC, theta0 = set_country_tau(best$theta, dC, log(7)), cores = cores, verbose = FALSE)
  c5 = jm_fit(dC, theta0 = set_country_tau(jm_theta0(dC, set_c), dC, log(3)), cores = cores, verbose = FALSE)
  c6 = jm_fit(dC, theta0 = set_country_tau(jm_theta0(dC, set_c), dC, log(14)), cores = cores, verbose = FALSE)
  best_of(a, b, c3, c4, c5, c6) })
f17 = cached("fit_tau_shared_R0_1.7.rds", {
  s = modifyList(jm_settings(), list(R0_fixed = 1.7)); d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  jm_fit(d, cores = cores, verbose = FALSE) })
# the pin reaches tau? tau by country at R0 = 1.7: every S0 should scale by 1.5/1.7 (the exact R0 x S0
# exchange) except where the ceiling released it, and the taus should not move
fC17 = cached("fit_tau_country_R0_1.7.rds", {
  s = modifyList(jm_settings(), list(R0_fixed = 1.7, tau_by_country = TRUE))
  d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  jm_fit(d, theta0 = jm_theta0(d, s), cores = cores, verbose = FALSE) })
f17_0 = cached("fit_tau0_R0_1.7.rds", {
  s = modifyList(jm_settings(), list(R0_fixed = 1.7, pr_tau_mean = -30, pr_tau_sd = 0.01))
  d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  jm_fit(d, theta0 = jm_theta0(d, s), cores = cores, verbose = FALSE) })

# ---- 2. profile over the shared tau ----
prof = cached("profile_tau.rds", {
  rows = data.frame(tau = c(0, exp(fS$theta[["log_tau"]])), loglik = c(f0$loglik, fS$loglik))
  for (tv in c(1, 4, 7, 10, 14, 21)){
    s = modifyList(jm_settings(), list(pr_tau_mean = log(tv), pr_tau_sd = 0.005))
    d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
    th0 = fS$theta; th0[["log_tau"]] = log(tv)
    f = jm_fit(d, theta0 = th0, cores = cores, verbose = FALSE)
    # the tau the fit actually held: the needle prior bends a little where the likelihood pulls hard
    rows = rbind(rows, data.frame(tau = exp(f$theta[["log_tau"]]), loglik = f$loglik))
  }
  list(profile = rows[order(rows$tau), ]) })

# ---- 3. per-country profile of tau_c, with that country's S0 along it ----
prof_c = cached("profile_tau_country.rds", {
  d = fC$d; th = fC$theta; grid = c(0, 1, 2, 4, 7, 10, 14, 21, 30, 45)
  do.call(rbind, parallel::mclapply(seq_len(d$n_country), function(ic){
    idx = d$off_country[ic] + seq_len(d$n_local[ic]); jt = max(idx); free = setdiff(idx, jt)
    mine = d$cs_of_country[[ic]] + 1L
    do.call(rbind, lapply(grid, function(tv){
      t2 = th; t2[jt] = if (tv == 0) -30 else log(tv)
      f = function(x){ t3 = t2; t3[free] = x; jm_country_negll_cpp(t3, d, ic - 1L) }
      o = optim(t2[free], f, method = "BFGS", control = list(maxit = 600, reltol = 1e-11)); t2[free] = o$par
      data.frame(country = d$countries[ic], tau = tv,
                 loglik = sum(jm_wave_loglik(list(theta = t2, d = d))$loglik[mine]),
                 S0 = plogis(t2[d$off_country[ic] + 1]))
    }))
  }, mc.cores = cores)) })

# ---- 3b. the theory plane: the model's own noise-free wave, S0 swept at fixed spreads ----
# rise = log(0.6/0.1) over the time the curve takes from 10% to 60% of its peak on the way up, width =
# full width at half maximum; both from interpolated crossing times on the WEEKLY curve, so they are
# continuous in S0 and measure what a weekly series can show
shape_theory = cached("shape_stats.rds", {
  th = f0$theta; d = f0$d; S = d$n_season
  ic = which(d$countries == "DK"); i_ref = d$cs_of_country[[ic]][1] + 1L; base = d$off_country[ic]
  xs = c(th[seq_len(S - 1)], -sum(th[seq_len(S - 1)])); s_ref = d$cs_season[i_ref] + 1L
  jI0 = base + 5 + d$n_src[ic] + d$cs_pos[i_ref] + 1L
  cross = function(v, lev, side){ p = which.max(v); H = v[p]
    k = if (side < 0) max(c(0, which(seq_len(length(v)) < p & v < lev * H))) else min(c(length(v) + 1, which(seq_len(length(v)) > p & v < lev * H)))
    if (k < 1 || k > length(v)) return(NA_real_)
    if (side < 0) k + (lev * H - v[k]) / (v[k + 1] - v[k]) else k - (lev * H - v[k]) / (v[k - 1] - v[k]) }
  shape = function(v) c(rise = log(6) / (cross(v, 0.6, -1) - cross(v, 0.1, -1)),
                        fwhm = cross(v, 0.5, 1) - cross(v, 0.5, -1))
  one = function(S0, tau){
    t2 = th; t2[base + 1] = qlogis(S0) - xs[s_ref]; t2[["log_tau"]] = if (tau == 0) -30 else log(tau)
    t2[base + 5 + seq_len(d$n_src[ic])] = -30                  # no off-season floor: the pure wave
    r = d$gamma * (d$R0_fixed * S0 - 1); t2[jI0] = log(1e-7)
    for (it in 1:4){ m = rowSums(jm_fitted_cpp(t2, d)$mu[[i_ref]])   # seed so the peak sits in week 25
      t2[jI0] = min(t2[jI0] + (which.max(m) - 25) * 7 * r, log(0.01)) }
    sh = shape(rowSums(jm_fitted_cpp(t2, d)$mu[[i_ref]]))
    data.frame(tau = tau, S0 = S0, rise = sh[["rise"]], fwhm = sh[["fwhm"]]) }
  list(theory = do.call(rbind, lapply(c(0, 7, 14, 21), function(tau)
    do.call(rbind, lapply(seq(0.74, 0.995, by = 0.005), function(S0) one(S0, tau)))))) })

# ---- 4. Hessians ----
idS = cached("ident_tau_shared.rds", jm_identifiability(fS, verbose = FALSE))
idC = cached("ident_tau_country.rds", jm_identifiability(fC, verbose = FALSE))

# ---- 5. recovery ----
# two more replicates per spread, so one noisy draw cannot pass for a bias
rec_more = cached("recovery_tau_more.rds", {
  rows = list()
  for (tv in c(7, 14)) for (rp in 1:2){
    truth = fS$theta; truth[["log_tau"]] = log(tv)
    r = jm_recover_once(dS, truth, seed = 4100L + 10L * round(tv) + rp, fit_args = list(cores = cores), with_intervals = FALSE)
    cm = r$comparison; s0 = cm[grepl(":logit_S0", cm$parameter), ]
    rows[[length(rows) + 1]] = data.frame(tau_true = tv, rep = rp, tau_est = exp(cm$est_theta[cm$parameter == "log_tau"]),
      S0_logit_bias = mean(s0$est_theta - s0$truth_theta), S0_rank = cor(s0$truth, s0$estimate, method = "spearman"),
      loglik = r$loglik)
  }
  do.call(rbind, rows) })
# is the low recovery the likelihood's own ridge or the priors? One 7-day data set, three fits: tau held
# at the truth, tau held at 3 days, and tau free with a threefold weaker S0 prior
ridge = cached("ridge_diagnosis.rds", {
  truth = fS$theta; truth[["log_tau"]] = log(7)
  ds = jm_simulate(truth, dS, seed = 3107L)
  one = function(pr_mean, pr_sd, S0_sd = 1){ dd = ds; dd$pr_tau_mean = pr_mean; dd$pr_tau_sd = pr_sd; dd$pr_S0_sd = S0_sd
    th0 = truth; th0[["log_tau"]] = if (pr_sd < 0.1) pr_mean else log(7)
    f = jm_fit(dd, theta0 = th0, cores = cores, verbose = FALSE)
    data.frame(tau = exp(f$theta[["log_tau"]]), loglik = f$loglik,
               S0_logit_error = mean(qlogis(vapply(jm_unpack(f$theta, dd)$country, `[[`, 0, "S0")) -
                                     qlogis(vapply(jm_unpack(truth, dS)$country, `[[`, 0, "S0")))) }
  cbind(fit = c("tau held at the true 7 d", "tau held at 3 d", "tau free, S0 prior sd 3"),
        rbind(one(log(7), 0.005), one(log(3), 0.005), one(log(7), 0.7, S0_sd = 3))) })
rec = cached("recovery_tau.rds", {
  res = list()
  for (tv in c(exp(fS$theta[["log_tau"]]), 7, 14)){
    truth = fS$theta; truth[["log_tau"]] = log(tv)
    res[[sprintf("shared_%g", tv)]] = list(rec = jm_recover_once(dS, truth, seed = 3100L + round(tv),
                                           fit_args = list(cores = cores), with_intervals = FALSE))
  }
  big = vapply(dC$N, sum, 0) > median(vapply(dC$N, sum, 0))
  truthC = fC$theta
  for (ic in seq_len(dC$n_country)) truthC[dC$off_country[ic] + dC$n_local[ic]] = log(if (big[ic]) 14 else 3)
  dsB = jm_simulate(truthC, dC, seed = 3201L)
  fB_c = jm_fit(dsB, theta0 = jm_theta0(dsB, set_c), cores = cores, verbose = FALSE)
  dsB_s = dS; dsB_s$y = dsB$y; dsB_s$lgamma_y1 = dsB$lgamma_y1; dsB_s$rates = dsB$rates
  fB_s = jm_fit(dsB_s, cores = cores, verbose = FALSE)
  res$country_hetero = list(fit_c = fB_c, fit_s = fB_s, truth = truthC, big = big,
                            compare = jm_compare_waves(fB_s, fB_c, "by-country vs shared, heterogeneous truth"))
  res })

# ---- 6. comparisons ----
cmp = cached("compare_tau.rds", rbind(
  jm_compare_waves(f0, fS, "shared tau vs none"),
  jm_compare_waves(fS, fC, "tau by country vs shared tau"),
  jm_compare_waves(f0, fC, "tau by country vs none")))
cat("=== fits ===\n")
for (nm in c("f0", "fS", "fC", "f17_0", "f17", "fC17")){ f = get(nm)
  cat(sprintf("  %-6s loglik %.2f  params %d  converged %s\n", nm, f$loglik, length(f$theta), f$converged)) }
cat("\n=== tau by country at R0 = 1.5 and 1.7 ===\n")
print(round(rbind(tau_1.5 = taus_of(fC), tau_1.7 = taus_of(fC17),
                  S0_1.5 = vapply(jm_unpack(fC$theta, fC$d)$country, `[[`, 0, "S0"),
                  S0_1.7 = vapply(jm_unpack(fC17$theta, fC17$d)$country, `[[`, 0, "S0")), 3))
cat("\n=== the shared-tau profile ===\n"); print(prof$profile, row.names = FALSE, digits = 7)
cat("\n=== recovery of a shared tau (extra replicates) and the ridge diagnosis ===\n")
print(rec_more, row.names = FALSE, digits = 3); print(ridge, row.names = FALSE, digits = 6)
cat("\n=== wave-level comparisons ===\n")
print(cmp[, c("contrast", "favour_B", "favour_A", "sign_p", "wilcoxon_p", "total", "total_deflated",
              "ci_wave_lo", "ci_wave_hi", "ci_country_lo", "ci_country_hi")], row.names = FALSE, digits = 3)
source("code/07_joint_model/plot_tau_analysis.R")
ggsave("output/joint_model/18_fatness.png",
       plot_tau_analysis(f0 = f0, fS = fS, fC = fC, prof = prof$profile, prof_c = prof_c, rec = rec,
                         rec_more = rec_more, cmp = cmp, theory = shape_theory$theory),
       width = 17, height = 13, dpi = 110)
cat("wrote output/joint_model/18_fatness.png\n")
