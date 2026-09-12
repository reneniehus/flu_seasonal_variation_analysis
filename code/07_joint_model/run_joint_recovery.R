# run_joint_recovery.R -- THE RECOVERY TEST. Simulate from a known truth onto the real design, refit
# the whole pipeline from scratch, and ask what came back.
#
# Three studies, in increasing ambition:
#   1. LOCAL RECOVERY      truth = the fitted optimum. Is the optimum recoverable, and do the 95%
#                          intervals actually cover at 95%?
#   2. PRIOR-SPACE RECOVERY truth drawn from the priors. Does recovery hold across the plausible
#                          parameter space, not just near one point?
#   3. DRIVER RECOVERY     truth in which a covariate really does move the season parameters by a
#                          known amount. Runs the two-step procedure the learning layer will use and
#                          checks the driver effect comes back. This is the learning layer's own
#                          validity test and the reason the harness is built to be pluggable.
#
#   Rscript code/07_joint_model/run_joint_recovery.R [n_local] [n_prior] [n_driver]
# Defaults 8 / 4 / 6. Each replicate is a full refit, so budget roughly two minutes per replicate
# plus four and a half for the Hessian wherever intervals are requested.
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_recovery.R")
source("code/07_joint_model/joint_report.R")
jm_load_cpp()

args = commandArgs(TRUE)
n_local  = if (length(args) >= 1) as.integer(args[1]) else 8L
n_prior  = if (length(args) >= 2) as.integer(args[2]) else 4L
n_driver = if (length(args) >= 3) as.integer(args[3]) else 6L
cores = max(1L, parallel::detectCores() - 1L)

stopifnot(file.exists("output/joint_model/joint_fit.rds"))
o = readRDS("output/joint_model/joint_fit.rds"); fit0 = o$fit; d = fit0$d
cat(sprintf("base fit: %d countries, %d country-seasons, %d parameters, negll %.1f\n",
            d$n_country, d$n_cs, d$n_par, fit0$negll))

out = list(base = list(negll = fit0$negll, loglik = fit0$loglik, seconds = fit0$seconds))

# ---- 1. local recovery, with intervals so coverage can be checked ----
if (n_local > 0){
  cat(sprintf("\n=== 1. local recovery: truth = the fitted optimum, %d replicates (with intervals) ===\n", n_local))
  rec1 = jm_recovery(d, truth_fn = function(i) jm_truth_from_fit(fit0), n_rep = n_local,
                     seed0 = 2000L, fit_args = list(cores = cores))
  s1 = jm_recovery_summary(rec1, d)
  cat("\n-- by parameter family --\n"); print(s1$by_family, row.names = FALSE, digits = 3)
  cat("\n-- rank recovery of the shared season parameters --\n"); print(s1$rank_shared, row.names = FALSE, digits = 3)
  cat(sprintf("\nsusceptibility RANKING across countries, median Spearman: %.3f\n", s1$rank_S0_spearman))
  cat(sprintf("median 95%% interval coverage across families: %.1f%%\n",
              100 * median(s1$by_family$coverage, na.rm = TRUE)))
  out$local = list(rec = rec1, summary = s1)
}

# ---- 2. prior-space recovery ----
if (n_prior > 0){
  cat(sprintf("\n=== 2. prior-space recovery: truth drawn from the priors, %d replicates ===\n", n_prior))
  rec2 = jm_recovery(d, truth_fn = function(i) jm_truth_from_prior(d, anchor = fit0$theta),
                     n_rep = n_prior, seed0 = 3000L,
                     fit_args = list(cores = cores))
  s2 = jm_recovery_summary(rec2, d)
  cat("\n-- by parameter family --\n"); print(s2$by_family, row.names = FALSE, digits = 3)
  cat("\n-- rank recovery of the shared season parameters --\n"); print(s2$rank_shared, row.names = FALSE, digits = 3)
  out$prior = list(rec = rec2, summary = s2)
}

# ---- 3. driver recovery: the learning layer's own validity test ----
if (n_driver > 0){
  cat(sprintf("\n=== 3. driver recovery: a covariate really moves the season parameters, %d replicates ===\n", n_driver))
  # a stand-in covariate with the shape a real driver would have (e.g. an H3N2-dominant indicator).
  # The point is not this particular x but whether a KNOWN effect of a covariate on the season
  # parameters survives the fit-then-regress procedure.
  x = c(1, 0, 1, 0, 0, 1, 0, 1)[seq_len(d$n_season)]
  beta_R0 = 0.06; beta_delta = 0.35
  cat(sprintf("truth: a one-sd move in the covariate shifts log R0 by %.2f (about %.0f%% on R0) and\n",
              beta_R0, 100 * (exp(beta_R0) - 1)))
  cat(sprintf("       log season visibility by %.2f (about %.0f%%)\n", beta_delta, 100 * (exp(beta_delta) - 1)))
  rec3 = jm_recovery(d, truth_fn = function(i) jm_truth_with_driver(d, x, beta_R0, beta_delta, anchor = fit0$theta),
                     n_rep = n_driver, seed0 = 4000L,
                     fit_args = list(cores = cores))
  s3 = jm_recovery_summary(rec3, d)
  dr = jm_driver_recovery(rec3, d)
  cat("\n-- did the driver effect come back? --\n"); print(dr$summary, row.names = FALSE, digits = 3)
  cat("\n-- by parameter family --\n"); print(s3$by_family, row.names = FALSE, digits = 3)
  out$driver = list(rec = rec3, summary = s3, driver = dr, x = x,
                    beta_R0 = beta_R0, beta_delta = beta_delta)
}

saveRDS(out, "output/joint_model/joint_recovery.rds")
cat("\n--- figure ---\n")
if (!is.null(out$local)) {
  ggplot2::ggsave("output/joint_model/13_recovery.png", plot_jm_recovery(out$local$rec, d, out$local$summary),
                  width = 10, height = 7, dpi = 115)
  cat("wrote output/joint_model/13_recovery.png\n")
}
cat("DONE\n")
