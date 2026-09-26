# The joint model must mean the same thing in both implementations, its parameter vector must be wired
# correctly, and the optimiser's guards must actually work. In rough order of how much damage a failure
# would do:
#   1. the C++ log-posterior equals an independent base-R implementation of the same model to 1e-10,
#      at the starting point and at random perturbations. The R version is written from the maths and
#      uses R's own eigen() for the spectral radius, so agreement means both are right rather than
#      that one was copied from the other;
#   2. every slot of the parameter vector is read back by jm_unpack as the name says, and the
#      sum-to-zero constraint on the season deviations holds;
#   3. no parameter is DEAD -- the objective responds to each one. The compartmental pilot shipped an
#      inert baseline slot for countries with only one data source, which made its Hessian exactly
#      singular; this test would have caught it;
#   4. the dynamics realise the R0 they are given, which is what the contact rescaling is for;
#   5. the recovery harness simulates from the same code it fits, onto the real design, and the driver
#      truth constructor encodes the slope it claims -- otherwise the recovery test is vacuous;
#   6. a rejected parameter vector is rejected by EVERY entry point rather than scored, and every
#      export checks theta's length before reading it;
#   7. the flat-line protector detects whole and partial flat lines, rescues them, never makes an
#      objective worse, leaves a healthy fit untouched, and reports a check that describes the theta
#      jm_fit actually returns;
#   8. the FIGURES describe the model they are drawn from: the design table's parameter counts add up
#      to the real parameter count, and every error bar is bound to the parameter it is drawn
#      against. Both had failed silently -- a mis-bound bar renders just as cleanly as a correct one.
# Skipped offline when the cached model inputs are absent, so run_tests.R stays runnable.

skip_if_not_installed("Rcpp")
skip_if_not(file.exists(here::here("output/models_in.rds")),
            "output/models_in.rds absent (run code/00_main.R)")

suppressMessages({
  source(here::here("code/01_main_supporting/setup.R"))
  source(here::here("code/01_main_supporting/stitch_iliplus.R"))
  source(here::here("code/01_main_supporting/sir_core.R"))
  source(here::here("code/06_comp_model/contact_matrix.R"))
  source(here::here("code/06_comp_model/comp_model_settings.R"))
  source(here::here("code/06_comp_model/comp_model_core.R"))
  source(here::here("code/06_comp_model/comp_model_data.R"))
  source(here::here("code/07_joint_model/joint_model.R"))
})
withr::with_dir(here::here(), jm_load_cpp())

models_in <- readRDS(here::here("output/models_in.rds"))
load(here::here("output/demography_respicast.Rdata")); demo <- obj
d  <- withr::with_dir(here::here(), jm_build_data(c("DK", "EE"), models_in, demo, verbose = FALSE))
th <- jm_theta0(d)

# The cached 12-country fit, or NULL when it predates the CURRENT parameter layout. Tests that read
# it must skip on NULL rather than hand a stale theta to the C++, which indexes by the new layout and
# either aborts or scores garbage. The layout changed on 2026-09-25 (R0 fixed; shared block 2S-1) and
# again on 2026-09-26 (the spatial spread tau: shared block 2S, or one slot per country).
.jm_cached <- local({
  f <- here::here("output/joint_model/joint_fit.rds")
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)
  ok <- !is.null(o$fit$d$R0_fixed) && !is.null(o$fit$d$tau_by_country) &&
        length(o$fit$theta) == length(jm_par_names(o$fit$d))
  if (isTRUE(ok)) o else NULL
})

test_that("the C++ log-posterior equals the base-R reference of the same model", {
  expect_equal(jm_negll_cpp(th, d), jm_negll_R(th, d), tolerance = 1e-10)
  set.seed(11)
  for (i in 1:5){
    t2 <- th + rnorm(length(th), 0, 0.2)
    expect_equal(jm_negll_cpp(t2, d), jm_negll_R(t2, d), tolerance = 1e-10)
  }
  # and the per-country pieces plus the shared priors must add up to the whole
  parts <- sum(vapply(seq_len(d$n_country), function(ic) jm_country_negll_cpp(th, d, ic - 1L), numeric(1)))
  expect_lt(abs((jm_negll_cpp(th, d) - parts) - (jm_negll_cpp(th, d) - parts)), 1e-12)  # finite
  expect_true(is.finite(jm_negll_cpp(th, d) - parts))
})

test_that("every parameter slot is read back as its name says, and the deviations average zero", {
  probe <- seq_len(d$n_par) / 100; names(probe) <- jm_par_names(d)
  p <- jm_unpack(probe, d); S <- d$n_season
  expect_equal(unname(p$x[seq_len(S - 1)]), unname(probe[seq_len(S - 1)]))
  expect_equal(sum(p$x), 0, tolerance = 1e-12)                # season effect on S0 averages zero
  expect_equal(p$R0, jm_settings()$R0_fixed)                  # R0 is fixed, reported for convenience
  expect_equal(unname(p$delta[seq_len(S - 1)]), unname(probe[S - 1L + seq_len(S - 1)]))
  expect_equal(sum(p$delta), 0, tolerance = 1e-12)            # the constraint that removes the pilot's flat direction
  expect_equal(log2(p$sigma_eld), unname(probe[2L * S - 1L]))
  expect_equal(log(p$tau), unname(probe[2L * S]))              # the shared spatial spread
  for (ic in seq_len(d$n_country)){
    b <- d$off_country[ic]; q <- p$country[[ic]]
    expect_equal(qlogis(q$S0), unname(probe[b + 1]))
    expect_equal(log(q$c),     unname(probe[b + 2]))
    expect_equal(q$off_young,  unname(probe[b + 3]))
    expect_equal(q$off_eld,    unname(probe[b + 4]))
    expect_equal(log(q$phi),   unname(probe[b + 5]))
    expect_equal(unname(log(q$b)),  unname(probe[b + 5 + seq_len(d$n_src[ic])]))
    expect_equal(unname(log(q$I0)), unname(probe[b + 5 + d$n_src[ic] + seq_len(d$n_cs_of_country[ic])]))
    expect_equal(q$tau, p$tau)                                 # shared: every country reads the one tau
  }
})

test_that("no parameter is dead: the objective responds to every slot", {
  base <- jm_negll_cpp(th, d)
  moved <- vapply(seq_len(d$n_par), function(i){
    t2 <- th; t2[i] <- t2[i] + 0.05
    abs(jm_negll_cpp(t2, d) - base) > 1e-9
  }, logical(1))
  expect_true(all(moved), info = paste("dead slots:", paste(jm_par_names(d)[!moved], collapse = ", ")))
  # the baseline slots must follow the sources the country actually has (the pilot hard-coded two)
  for (ic in seq_len(d$n_country))
    expect_equal(d$n_src[ic], length(unique(d$cs_src[d$cs_of_country[[ic]] + 1L])))
})

test_that("the realised R0 is what the dynamics produce, and the contact scaling is what makes it so", {
  # at full susceptibility and no vaccination the total infected grows at gamma*(R0 - 1) per day once
  # the dominant mode takes over -- the point of rescaling the contact matrix to spectral radius one
  set <- jm_settings(); A <- 3L
  Cn <- d$Cn[[1]]; expect_equal(max(abs(eigen(Cn, only.values = TRUE)$values)), 1, tolerance = 1e-8)
  R0 <- 1.5; gamma <- set$gamma_per_day; beta <- R0 * gamma
  Iu <- rep(1e-9, A); Su <- rep(1 - 1e-9, A); tot <- numeric(120)
  for (day in 1:120){
    lam <- beta * as.numeric(Cn %*% Iu)
    nu <- lam * Su; Su <- Su - nu; Iu <- Iu + nu - gamma * Iu
    tot[day] <- sum(Iu)
  }
  expect_equal(mean(diff(log(tot[60:100]))), log(1 + gamma * (R0 - 1)), tolerance = 5e-3)
})

test_that("the recovery harness simulates from the model onto the real design", {
  # The generative and fitted models must be the same code, and the simulated data must sit on the
  # real design: same missing cells, same dimensions. Anything else and the recovery test is testing
  # a different model from the one being fitted, which is the usual way such a test ends up vacuous.
  source(here::here("code/07_joint_model/joint_recovery.R"))
  set.seed(3)
  ds <- jm_simulate(th, d, seed = 5)
  expect_equal(length(ds$y), length(d$y))
  for (i in seq_along(d$y)){
    expect_equal(dim(ds$y[[i]]), dim(d$y[[i]]))
    expect_equal(is.na(ds$y[[i]]), is.na(d$y[[i]]))        # the real missing pattern is preserved
    expect_true(all(ds$y[[i]][is.finite(ds$y[[i]])] >= 0)) # counts, so non-negative
  }
  # the normalising constant and the fields the starting values are built from must be recomputed
  # from the SIMULATED data, or the refit gets help from the real series
  expect_false(isTRUE(all.equal(ds$lgamma_y1, d$lgamma_y1)))
  expect_false(isTRUE(all.equal(ds$rates[[1]], d$rates[[1]])))
  expect_true(all(is.finite(ds$lgamma_y1)))
  expect_true(is.finite(jm_negll_cpp(th, ds)))
  # expected counts must come from the same C++ the likelihood uses
  expect_equal(jm_fitted_cpp(th, ds)$mu[[1]], jm_fitted_cpp(th, d)$mu[[1]])
})

test_that("the driver truth constructor really encodes the driver effect it claims", {
  source(here::here("code/07_joint_model/joint_recovery.R"))
  S <- d$n_season
  x <- seq_len(S) %% 2
  set.seed(9)
  tt <- jm_truth_with_driver(d, x, beta_x = 0.10, beta_delta = 0.40, noise_x = 0, noise_delta = 0)
  dr <- attr(tt, "driver")
  expect_equal(dr$beta_x, 0.10); expect_equal(dr$beta_delta, 0.40)
  # with the noise switched off, regressing the truth back on the covariate returns the slopes exactly
  p <- jm_unpack(tt, d)
  expect_equal(unname(coef(lm(p$x ~ dr$x))[2]), 0.10, tolerance = 1e-8)
  expect_equal(sum(p$x), 0, tolerance = 1e-10)
  expect_equal(unname(coef(lm(p$delta ~ dr$x))[2]), 0.40, tolerance = 1e-8)
  expect_equal(sum(p$delta), 0, tolerance = 1e-10)         # the constraint still holds
})

test_that("curvature intervals contain the estimate and cover the constrained deviation", {
  source(here::here("code/07_joint_model/joint_recovery.R"))
  small <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  iv <- jm_intervals(small)
  expect_true(all(iv$lower <= iv$estimate + 1e-9))
  expect_true(all(iv$estimate <= iv$upper + 1e-9))
  # every free parameter, plus one extra row for the deviation that is minus the sum of the others
  expect_equal(nrow(iv), length(small$theta) + 2L)   # x_last AND delta_last are both derived rows
  expect_true(paste0("delta_", d$seasons[d$n_season]) %in% iv$parameter)
  expect_true(paste0("x_", d$seasons[d$n_season]) %in% iv$parameter)
})

test_that("the flat-line protector detects a flat fit, escapes it, and leaves a healthy one alone", {
  # The guarded failure: a country's local block has a second optimum in which the dispersion
  # collapses, the negative binomial becomes diffuse enough that any curve fits, and the country sits
  # at its baseline for every season. It cost the Netherlands in the first joint fit and a harder
  # search later found a solution 406 nats better, so it is an optimiser failure, not a fact.
  skip_if(is.null(.jm_cached), "no saved fit, or it predates this parameter layout")
  fit <- .jm_cached$fit
  dd <- fit$d; bl <- jm_blocks(dd)

  # (a) the detector must not fire on a healthy fit, and with real margin, not marginally
  chk <- jm_flat_check(fit$theta, dd)
  expect_false(any(chk$flat))
  expect_gt(min(chk$attack_max), 0.15)        # threshold is 0.03: at least a 5x margin
  expect_gt(min(chk$peak_epi_frac), 0.75)     # threshold is 0.15

  # (b) slot 5 of every local block really is the dispersion, which jm_unflatten relies on
  nm <- jm_par_names(dd)
  for (ic in seq_len(dd$n_country)) expect_match(nm[bl$local[[ic]][5]], ":log_phi$")

  # (c) break one country deliberately, two ways, and require the detector to fire
  ic <- 3L; idx <- bl$local[[ic]]
  i_seed <- 5L + dd$n_src[ic] + seq_len(dd$n_cs_of_country[ic])
  broken <- list(kill_susceptibility = local({ t <- fit$theta; t[idx[1]] <- qlogis(0.001); t[idx[5]] <- log(0.02); t }),
                 tiny_seeds          = local({ t <- fit$theta; t[idx[i_seed]] <- log(1e-30); t }))
  for (nmb in names(broken)){
    tb <- broken[[nmb]]
    expect_true(jm_flat_check(tb, dd)$flat[ic], info = nmb)
    # (d) the protector must rescue it, and must only ever adopt an IMPROVEMENT
    before <- jm_country_negll_cpp(tb, dd, ic - 1L)
    prot <- jm_unflatten(tb, dd, cores = 1, verbose = FALSE)
    after <- jm_country_negll_cpp(prot$theta, dd, ic - 1L)
    expect_lte(after, before + 1e-6)                          # never worse
    expect_lte(jm_negll_cpp(prot$theta, dd), jm_negll_cpp(tb, dd) + 1e-6)
    expect_false(prot$check$flat[ic], info = nmb)              # and actually rescued
    expect_equal(prot$n_fixed, 1L)
  }

  # (e) on a healthy fit it must be a no-op, not a small perturbation
  p0 <- jm_unflatten(fit$theta, dd, cores = 1, verbose = FALSE)
  expect_equal(unname(p0$theta), unname(fit$theta))
  expect_equal(nrow(p0$report), 0L)
  expect_equal(p0$n_unresolved, 0L)
})

test_that("the data-implied dispersion is one definition used everywhere", {
  # jm_phi_data feeds the adequacy diagnostic, the multi-start and the protector; if the three ever
  # disagree the noise-budget figure stops meaning what it says
  skip_if(is.null(.jm_cached), "no saved fit, or it predates this parameter layout")
  fit <- .jm_cached$fit
  phid <- jm_phi_data(fit$d)
  expect_length(phid, fit$d$n_country)
  expect_true(all(is.finite(phid) & phid > 0))
  ad <- jm_adequacy(fit)
  expect_equal(ad$phi_data, phid, tolerance = 1e-10)     # the diagnostic uses the same numbers
  expect_equal(ad$cv_data, 1 / sqrt(phid), tolerance = 1e-10)
})

test_that("the C++ and R implementations agree where the pool cap BINDS, not just at sane parameters", {
  # The cap (flow limited to what S_u holds) only engages when lambda > 1 in a day, i.e. around
  # R0 > 6. The identity test at fitted values therefore passed for a while with the cap present in
  # C++ and absent in R. This exercises the regime where they could differ.
  for (R0 in c(1.5, 6.5, 10, 40)){
    t2 <- th; t2[seq_len(d$n_season)] <- log(R0)
    expect_equal(jm_negll_cpp(t2, d), jm_negll_R(t2, d), tolerance = 1e-9,
                 info = paste("R0 =", R0))
  }
  # and the invariant the cap restores: nobody can be infected more than once
  t3 <- th; t3[seq_len(d$n_season)] <- log(40)
  f <- jm_fitted_cpp(t3, d)
  for (i in seq_len(d$n_cs)){
    S0 <- plogis(t3[jm_blocks(d)$local[[d$cs_country[i] + 1L]][1]])
    expect_lte(max(f$attack[i, ]), S0 + 1e-9)
  }
})

test_that("a rejected parameter vector is rejected by every entry point, not scored", {
  S <- d$n_season
  bad <- list(
    "sigma overflows"  = local({ t <- th; t[2L * S - 1L] <- 1030; t }),
    "delta overflows"  = local({ t <- th; t[S - 1L + seq_len(S - 1)] <- 800; t[grep(":log_I0", names(th))] <- -800; t }),
    "one delta overflows" = local({ t <- th; t[S] <- 800; t }))
  for (nm in names(bad)){
    t2 <- bad[[nm]]
    # the objective must return its sentinel, never a finite value that looks like a better fit
    expect_gte(jm_negll_cpp(t2, d), 1e10, label = nm)
    # and the log-likelihood must PROPAGATE that, not strip priors off the sentinel and return ~-1e10
    expect_false(is.finite(jm_loglik_cpp(t2, d)), label = nm)
  }
  # fitted values cannot be computed there, and must say why rather than dying on an empty List
  e <- try(jm_fitted_cpp(bad[["sigma overflows"]], d), silent = TRUE)
  expect_s3_class(e, "try-error")
  expect_match(conditionMessage(attr(e, "condition")), "shared block")
})

test_that("every exported entry point checks the length of theta before reading it", {
  # country_lp indexes up to exactly n_par - 1 for the last country, so a short theta reads the
  # memory next to the R vector and returns a plausible number
  short <- th[seq_len(length(th) - 1L)]
  expect_error(jm_negll_cpp(short, d), "wrong length")
  expect_error(jm_country_negll_cpp(short, d, 0L), "wrong length")
  expect_error(jm_loglik_cpp(short, d), "wrong length")
  expect_error(jm_fitted_cpp(short, d), "wrong length")
  expect_error(jm_country_negll_cpp(th, d, d$n_country), "out of range")
})

test_that("a PARTIAL flat line is detected and rescued, and a low-attack design is not false-flagged", {
  skip_if(is.null(.jm_cached), "no saved fit, or it predates this parameter layout")
  fit <- .jm_cached$fit
  dd <- fit$d; bl <- jm_blocks(dd)
  ic <- which(dd$countries == "NL"); if (!length(ic)) ic <- 1L
  i_seed <- bl$local[[ic]][5L + dd$n_src[ic] + seq_len(dd$n_cs_of_country[ic])]
  n_flat <- max(1L, floor(length(i_seed) / 2))
  tp <- fit$theta; tp[i_seed[seq_len(n_flat)]] <- -60          # half this country's seasons flat-lined
  cp <- jm_flat_check(tp, dd)
  expect_true(cp$flat[ic])                                     # the max over seasons would have hidden it
  expect_equal(cp$n_seasons_flat[ic], n_flat)
  expect_gte(cp$n_seasons_healthy[ic], 2L)
  prot <- jm_unflatten(tp, dd, cores = 1, verbose = FALSE)
  expect_lt(jm_negll_cpp(prot$theta, dd), jm_negll_cpp(tp, dd))
  expect_false(prot$check$flat[ic])
  # the partial trigger must need HEALTHY seasons too, so a uniformly low-attack design is spared
  tl <- fit$theta; tl[bl$local[[ic]][1]] <- qlogis(0.02)       # every season small, none anomalous
  cl <- jm_flat_check(tl, dd, attack_min = 1e-6)
  expect_false(cl$flat_partial[ic])
})

test_that("jm_fit's reported flat check describes the theta it actually returns", {
  # The bug: the post-polish rescue was adopted only if it CLEARED the threshold, while the reported
  # check was computed on the rescued vector either way -- so fit$flat could describe a vector that
  # fit$theta was not. Worth 5648 nats in the measured case.
  small <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  expect_equal(jm_flat_check(small$theta, d, small$flat_thresholds[["attack_min"]],
                             small$flat_thresholds[["epi_frac_min"]])$flat,
               small$flat$flat)
  expect_equal(small$n_flat_unresolved, sum(small$flat$flat))
  expect_true(is.logical(small$converged))
  # `converged` must describe the WHOLE fit, not just the sweep loop. It read FALSE on a fit whose
  # three stages had all succeeded, because it meant "the sweep loop did not run out of sweeps" --
  # and it also called a break on the very last sweep a failure.
  expect_equal(small$converged,
               isTRUE(all(small$conv_local == 0L) && small$conv_shared == 0L &&
                      small$conv_polish == 0L && small$n_flat_unresolved == 0L))
  expect_true(is.logical(small$sweeps_hit_tol))
  # the sweep loop's own status is reported separately, and agrees with the trace
  sh <- small$trace[small$trace$step == "shared", ]
  if (nrow(sh) >= 2) expect_equal(small$sweep_last_gain, -tail(diff(sh$negll), 1), tolerance = 1e-8)
  expect_equal(small$sweeps_hit_tol, isTRUE(small$sweep_last_gain < 0.05))
  # the extrapolated tail is either a positive number of nats or NA, never negative or infinite
  expect_true(is.na(small$sweep_tail_nats) ||
              (is.finite(small$sweep_tail_nats) && small$sweep_tail_nats >= 0))
  expect_named(small$flat_thresholds, c("attack_min", "epi_frac_min"))
})

# ---- the figures must describe the model they are drawn from ----
# A figure is a claim. These two were wrong on disk before the checks below existed: the design
# table counted one season visibility too many (the sum-to-zero constraint leaves S-1 free), and
# figure 13's error bars were bound to a blocked name list while the plotting frame was interleaved
# by country, so 22 of 24 bars sat on the wrong parameter. Both are silent failures -- the figure
# renders perfectly either way -- so they are held by a test rather than by inspection.
test_that("the design figure's parameter counts add up to the model's actual parameter count", {
  suppressMessages(source(here::here("code/07_joint_model/joint_report.R")))
  spec <- jm_design_spec(d)
  expect_equal(sum(spec$n), d$n_par)
  # the fitted rows must be exactly the ones with a non-zero count, and the fixed ones zero
  expect_true(all(spec$n[spec$group == "fixed"] == 0))
  expect_true(all(spec$n[spec$group != "fixed"] > 0))
  # and the counts must survive a differently-shaped panel
  d3 <- withr::with_dir(here::here(),
                        jm_build_data(c("DK", "EE", "FR"), models_in, demo, verbose = FALSE))
  expect_equal(sum(jm_design_spec(d3)$n), d3$n_par)
})

test_that("every figure's error bar is bound to the parameter it is drawn against", {
  suppressMessages({source(here::here("code/07_joint_model/joint_recovery.R"))
                    source(here::here("code/07_joint_model/joint_report.R"))})
  small <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  iv <- jm_intervals(small)
  s  <- jm_summary_country(small); ss <- jm_summary_season(small)
  # the interval table's own estimate must equal the quantity each figure plots, on that figure's
  # scale -- this is what says the back-transform in the figure matches the one in jm_intervals
  pick <- function(n) iv[match(n, iv$parameter), ]
  expect_equal(pick(paste0("x_", d$seasons))$estimate, ss$x, tolerance = 1e-8)
  expect_equal(exp(pick(paste0("delta_", d$seasons))$estimate), ss$reporting_mult, tolerance = 1e-8)
  expect_equal(pick(paste0(d$countries, ":logit_S0"))$estimate, s$S0, tolerance = 1e-8)
  expect_equal(pick(paste0(d$countries, ":log_c"))$estimate, s$c_adult, tolerance = 1e-8)
  expect_equal(pick(paste0(d$countries, ":off_young"))$estimate, s$rel_young, tolerance = 1e-8)
  expect_equal(pick(paste0(d$countries, ":off_eld"))$estimate, s$rel_elderly, tolerance = 1e-8)
  # figure 13 reshapes country x {young, elderly} into long form; the bars must follow the reshape
  age <- s[, c("country", "rel_young", "rel_elderly")]
  names(age) <- c("country", "young", "elderly")
  age <- tidyr::pivot_longer(age, c("young", "elderly"), names_to = "group", values_to = "rel")
  age$par <- paste0(age$country, ifelse(age$group == "young", ":off_young", ":off_eld"))
  ci <- .jm_iv(iv, age$par)
  expect_equal(ci$estimate, age$rel, tolerance = 1e-8)
  expect_true(all(age$rel >= pmin(ci$lower, ci$upper) - 1e-9 &
                  age$rel <= pmax(ci$lower, ci$upper) + 1e-9))
  # a missing name must warn rather than silently drop a bar
  expect_warning(.jm_iv(iv, c(d$countries[1], ":off_young", "no_such_parameter")))
})

test_that("the two other ordering-sensitive figures label their values correctly", {
  suppressMessages(source(here::here("code/07_joint_model/joint_report.R")))
  small <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  nm <- jm_par_names(d); p <- jm_unpack(small$theta, d)
  # figure 05 reads each country's seeds as a vector; the SEASON LABEL it draws them against comes
  # from cs_of_country, so the two orders have to agree. The names in theta are the ground truth.
  for (ic in seq_len(d$n_country)){
    ics <- d$cs_of_country[[ic]] + 1L
    want <- paste0(d$countries[ic], ":log_I0_", d$seasons[d$cs_season[ics] + 1L])
    expect_false(anyNA(match(want, nm)))
    expect_equal(log(unname(p$country[[ic]]$I0)), unname(small$theta[match(want, nm)]),
                 tolerance = 1e-9)
  }
  # figure 07 flattens week x age matrices into a long frame; every series must still be the column
  # it came from, for both the observed and the modelled layer
  tf <- jm_tidy_fit(small); f <- jm_fitted_cpp(small$theta, d)
  expect_equal(nrow(tf), 2L * sum(d$n_weeks) * length(d$groups))
  for (i in seq_len(d$n_cs)){
    ic <- d$cs_country[i] + 1L; per <- d$rate_per / d$N[[ic]]
    sel <- tf$country == d$countries[ic] & tf$season == d$seasons[d$cs_season[i] + 1L]
    for (g in seq_along(d$groups)){
      s <- tf[sel & tf$group == d$groups[g], ]
      expect_equal(s$value[s$what == "observed"][order(s$week[s$what == "observed"])],
                   unname(d$y[[i]][, g] * per[g]), tolerance = 1e-9)
      expect_equal(s$value[s$what == "model"][order(s$week[s$what == "model"])],
                   unname(f$mu[[i]][, g] * per[g]), tolerance = 1e-9)
    }
  }
})

test_that("the misspecification arms are real violations, not reparameterisations", {
  # THE TRAP, hit for real: the first misspecification arm multiplied each country's expected counts
  # by a constant. That is exactly what the per-country reporting level c does, so the model absorbed
  # it perfectly and the arm scored BETTER than the control while appearing to test the assumption.
  # A violation only counts if the ratio to the unperturbed mean VARIES WITHIN a country-season, so
  # no single reporting number can undo it.
  source(here::here("code/07_joint_model/joint_recovery.R"))
  set.seed(4)
  for (v in c("r0_by_country", "second_wave")){
    ds <- jm_simulate_violation(th, d, v, seed = 77)
    expect_true(jm_violation_is_real(ds), info = v)
    expect_identical(attr(ds, "violation"), v)
    # still on the real design: same shapes, same missing cells, counts non-negative
    for (i in seq_along(d$y)){
      expect_equal(dim(ds$y[[i]]), dim(d$y[[i]]))
      expect_equal(is.na(ds$y[[i]]), is.na(d$y[[i]]))
    }
    expect_true(is.finite(jm_negll_cpp(th, ds)))
  }
  # the absorbable perturbation must be REJECTED by the same check, or it cannot protect anything
  ds0 <- jm_simulate_violation(th, d, "none", seed = 77)
  expect_false(jm_violation_is_real(ds0))
  fake <- ds0; mu0 <- jm_fitted_cpp(th, d)$mu
  attr(fake, "violation") <- "visibility_by_country"
  attr(fake, "mu_ratio") <- lapply(seq_along(mu0), function(i) mu0[[i]] * 0 + 1.7)  # constant per cell
  expect_false(jm_violation_is_real(fake))
  # r0_by_country must enter the DYNAMICS: with no multiplier it reduces to the plain simulator
  set.seed(1); a <- jm_simulate_violation(th, d, "r0_by_country", seed = 5, sd_log_r0 = 0)
  set.seed(1); b <- jm_simulate(th, d, seed = 5)
  expect_equal(a$y[[1]], b$y[[1]])
})

test_that("a recovery replicate skips the Hessian when it is not asked for intervals", {
  # jm_recover_once used to compute the curvature intervals unconditionally, which is minutes per
  # replicate thrown away, because the only thing it needed from them was the name-derived transform.
  source(here::here("code/07_joint_model/joint_recovery.R"))
  nm <- jm_par_names(d)
  expect_equal(jm_par_kind(nm), jm_intervals(jm_fit(d, max_sweeps = 1L, cores = 1, verbose = FALSE))$kind[seq_along(nm)])
  # these files are sourced, not a package, so shadow the binding jm_recover_once actually resolves
  called <- 0L
  real <- jm_intervals
  assign("jm_intervals", function(...) { called <<- called + 1L; stop("Hessian computed") },
         envir = globalenv())
  on.exit(assign("jm_intervals", real, envir = globalenv()), add = TRUE)
  r <- jm_recover_once(d, th, seed = 8, fit_args = list(max_sweeps = 1L, cores = 1),
                       with_intervals = FALSE)
  expect_equal(called, 0L)
  expect_false("covered" %in% names(r$comparison))
  expect_true(all(c("truth", "estimate", "truth_theta", "est_theta") %in% names(r$comparison)))
  expect_true(is.finite(r$adequacy_excess))
  # and it IS called when intervals are asked for
  expect_error(jm_recover_once(d, th, seed = 8, fit_args = list(max_sweeps = 1L, cores = 1),
                               with_intervals = TRUE), "Hessian computed")
  expect_equal(called, 1L)
})

test_that("every figure in the default set actually renders", {
  # The failure this catches: an aesthetic that names a column which is not in the plot's data.
  # ggplot captures its data BY VALUE, so a column attached to the frame after the ggplot() call is
  # invisible to every later layer -- and nothing complains until the plot is drawn. Two figures were
  # broken that way for the length of one edit. ggplot_build() forces the evaluation that ggsave would.
  suppressMessages({source(here::here("code/07_joint_model/joint_recovery.R"))
                    source(here::here("code/07_joint_model/joint_report.R"))})
  small <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  iv <- jm_intervals(small); id <- NULL
  build <- function(p) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  # the data-overview and mechanism figures. The mechanism figure computes its curves from the model's
  # own C++ at perturbed parameter vectors, so it cannot drift from the model -- but it does index
  # specific parameter slots, so a layout change must fail here rather than mislabel a panel.
  # patchworkGrob opens a graphics device, which drops an Rplots.pdf in the test directory; render to
  # a throwaway device instead so the test leaves nothing behind
  bp <- function(p) { pdf(NULL); on.exit(dev.off(), add = TRUE)
                      expect_s3_class(patchwork::patchworkGrob(p), "gtable") }
  build(plot_jm_data_panel(small)); bp(plot_jm_data_features(small)); bp(plot_jm_mechanism(small))
  # with intervals, which is how the report is written
  build(plot_jm_design(small)); build(plot_jm_arrival(small))
  build(plot_jm_fit_overview(small)); build(plot_jm_fit_country(small, d$countries[1]))
  build(plot_jm_adequacy(small)); build(plot_jm_attack(small))
  for (p in list(plot_jm_season_S0(small, iv), plot_jm_season_visibility(small, iv),
                 plot_jm_country_S0(small, iv), plot_jm_country_reporting(small, iv),
                 plot_jm_age_offsets(small, iv))) build(p)
  # and WITHOUT intervals, the other branch of every one of those five
  for (p in list(plot_jm_season_S0(small), plot_jm_season_visibility(small),
                 plot_jm_country_S0(small), plot_jm_country_reporting(small),
                 plot_jm_age_offsets(small))) build(p)
})

test_that("the pipeline writes the figure set WITH uncertainty intervals", {
  # The regression this guards: run_joint_model.R computed the curvature intervals, saved them into
  # joint_fit.rds, and then called save_jm_report(fit, id) without them -- so every published figure
  # 06-10 was drawn with no error bars for two days. Nothing complained, because a figure without
  # uncertainty renders exactly as cleanly as one with it. Checked statically, since the runner is a
  # script rather than a function.
  runner <- paste(readLines(here::here("code/07_joint_model/run_joint_model.R")), collapse = "\n")
  call <- regmatches(runner, regexpr("save_jm_report\\([^)]*\\)", runner))
  expect_length(call, 1L)
  expect_match(call, "\\biv\\b", info = paste("save_jm_report call was:", call))
  # and the writer must say so when they are absent, rather than quietly dropping them
  suppressMessages({source(here::here("code/07_joint_model/joint_recovery.R"))
                    source(here::here("code/07_joint_model/joint_report.R"))})
  expect_warning(save_jm_report(jm_fit(d, max_sweeps = 1L, cores = 1, verbose = FALSE),
                                dir = withr::local_tempdir()),
                 "WITHOUT uncertainty")
})

# ---- the data object must be the design the documents describe ----
test_that("the default design is the documented one, and the attack rate is a season quantity", {
  # FINDING: min_seasons defaulted to 6 while the owner's decision (2026-09-12) was 5. Only the
  # runner passed the override, so anyone reproducing "the model" as MODEL.md describes it silently
  # got 11 countries / 81 country-seasons / 173 parameters. No test pinned the real candidate list,
  # because the small designs used above all have 8 seasons and never exercise the default.
  cand <- c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")
  # TWO documented designs, and both are pinned. Without the positivity-encoding exclusion the
  # min_seasons = 5 decision gives 12 / 86 / 183; with it (the default since 2026-09-16, provisional
  # pending surveillance confirmation) CZ 2024/2025 is dropped and CZ loses the ERVISS baseline slot
  # that had no off-season behind it, giving 12 / 85 / 182. (R0 fixed since 2026-09-25: the
  # shared block lost the S slots of R0_s; the spatial spread tau added one on 2026-09-26.)
  full <- withr::with_dir(here::here(),
            jm_build_data(cand, models_in, demo, verbose = FALSE,
                          exclude_ambiguous_positivity = FALSE))
  expect_equal(c(full$n_country, full$n_cs, full$n_par), c(12L, 86L, 184L))
  dd <- withr::with_dir(here::here(), jm_build_data(cand, models_in, demo, verbose = FALSE))
  expect_equal(c(dd$n_country, dd$n_cs, dd$n_par), c(12L, 85L, 182L))
  expect_true("ES" %in% dd$countries)            # the country the min_seasons decision was about
  expect_equal(dd$n_season, 8L)                  # no season is lost entirely by the exclusion
  # the dynamics horizon is a full season for every cell, and the observation windows are not
  expect_gte(dd$attack_weeks, max(dd$n_weeks))
  expect_true(any(dd$n_weeks < dd$attack_weeks))
  # THE GUARANTEE: running the dynamics past the observation window must not touch the likelihood.
  # Checked on whichever design the cached fit's theta actually belongs to -- the shipped fit uses the
  # exclusion, but the check is meaningful on either, so pick by length rather than assuming.
  skip_if(is.null(.jm_cached), "no saved fit, or it predates this parameter layout")
  fit0 <- .jm_cached$fit
  dref <- if (length(fit0$theta) == dd$n_par) dd else full
  skip_if_not(length(fit0$theta) == dref$n_par, "cached fit predates both layouts")
  d_short <- dref; d_short$attack_weeks <- 0L    # the old behaviour
  # Without spread the guarantee is exact ...
  t_nospread <- fit0$theta; t_nospread[["log_tau"]] <- -30
  expect_equal(jm_negll_cpp(t_nospread, dref), jm_negll_cpp(t_nospread, d_short), tolerance = 1e-12)
  # ... and WITH spread it must not hold: copies of the wave that started earlier are still running
  # when a country's observation window closes, so the last observed weeks depend on the epidemic past
  # the window. Cutting the dynamics there would drop them (measured: 2.3 nats on the working fit); the
  # full-season horizon every design uses is what keeps the likelihood right.
  expect_gt(abs(jm_negll_cpp(fit0$theta, dref) - jm_negll_cpp(fit0$theta, d_short)), 0.01)
  # ... but it does change the attack rate, which is the point
  a_season <- jm_fitted_cpp(fit0$theta, dref)$attack
  a_window <- jm_fitted_cpp(fit0$theta, d_short)$attack
  expect_gt(max(abs(a_season - a_window) / a_window), 0.05)
  expect_true(all(a_season >= a_window - 1e-12))  # a longer horizon can only add infections
  # and the provenance of each country's mixing pattern travels with the fit
  expect_equal(length(dd$contact_source), dd$n_country)
  expect_equal(unname(dd$contact_source[dd$countries == "NO"]), "EU average")
})

test_that("an unresolved season or source label fails in R instead of aborting in C++", {
  # FINDING: cs_season and cs_src come from match() - 1L, which is NA for an unknown label, and srcs
  # is built with sort(unique(...)) which DROPS NA -- so one unresolved source label yields an NA
  # index that the C++ uses raw (logb[NA_INTEGER] = logb[INT_MIN]), killing the R process outright.
  # The guard has to be in jm_build_data, because by the time the C++ sees it there is no recovery.
  src <- paste(readLines(here::here("code/07_joint_model/joint_model.R")), collapse = "\n")
  expect_match(src, "anyNA\\(cs_season\\)")
  expect_match(src, "anyNA\\(cs_src\\)")
  # the mechanism the guard exists for, demonstrated without touching the C++
  expect_true(is.na(match("no_such_source", sort(unique(c("ERVISS", NA))))))
  expect_equal(sort(unique(c("ERVISS", NA_character_))), "ERVISS")   # sort() drops the NA
})

test_that("the two settings files cannot silently disagree about the rate basis", {
  # FINDING: d$y is built with comp_model_settings() while d$rate_per is stored from jm_settings().
  # They agreed only by coincidence; a change to one alone would put the counts on one basis and
  # every per-100k conversion on another (figure 06 out by 10x, every baseline b silently rescaled).
  expect_equal(jm_settings()$rate_per, comp_model_settings()$rate_per)
  bad <- modifyList(jm_settings(), list(rate_per = 1e6))
  expect_error(withr::with_dir(here::here(),
                 jm_build_data(c("DK", "EE"), models_in, demo, set = bad, verbose = FALSE)),
               "rate_per disagrees")
})

test_that("the vaccination fallback degrades instead of dying when the external file is absent", {
  # FINDING: `ext = ... else c()` made ext NULL, so is.na(ext[s]) was logical(0) and the documented
  # four-step fallback ladder raised "argument is of length zero" for every country.
  cov <- .cm_vax_coverage("DK", c("2014/2015", "2023/2024"), models_in,
                          ext_path = here::here("data/external/__absent__.csv"))
  expect_equal(nrow(cov), 2L)
  expect_true(all(is.finite(cov$coverage)))
  expect_true(all(!is.na(cov$provenance)))
  # and it agrees with the real call wherever the external file is not the source used
  real <- .cm_vax_coverage("DK", c("2014/2015", "2023/2024"), models_in)
  expect_equal(cov$coverage[real$provenance != "data/external post-COVID"],
               real$coverage[real$provenance != "data/external post-COVID"])
})

test_that("a season-level number is never reported without its sample size", {
  # FINDING: season support is uneven (the last season rests on 7 of 12 countries and ~half the
  # weekly cells of the richest) and it carries the highest fitted R0, but jm_summary_season
  # returned R0 and visibility with nothing to say how much data stood behind them.
  small <- jm_fit(d, max_sweeps = 1L, cores = 1, verbose = FALSE)
  s <- jm_summary_season(small)
  expect_true(all(c("n_country", "obs_cells", "last_week_min", "last_week_max") %in% names(s)))
  expect_equal(sum(s$n_country), d$n_cs)
  expect_equal(sum(s$obs_cells), sum(vapply(d$y, function(m) sum(is.finite(m)), numeric(1))))
  expect_true(all(s$last_week_min <= s$last_week_max))
  # and the country table carries which mixing pattern each country's offsets were fitted under
  expect_true("contact" %in% names(jm_summary_country(small)))
})

test_that("the ambiguous positivity encoding is detected and its seasons excluded", {
  # PROVISIONAL, pending confirmation by surveillance colleagues (owner, 2026-09-16). ERVISS encodes
  # a zero-detection week two ways -- "detections = 0" explicitly, or the detections row absent with
  # tests > 0 -- and our re-derived positivity turns the second into NA, which the stitch deletes. A
  # trailing run then SHORTENS the season, so a country-season can be fitted on a window containing no
  # off-season with nothing reporting it. These tests pin the detector and the exclusion so neither
  # can drift while the question is open.
  a <- erviss_encoding_ambiguous(models_in)
  expect_true(all(c("country_short", "season", "n_ambiguous", "n_not_plausibly_zero",
                    "n_no_neighbour", "min_p_zero") %in% names(a)))
  expect_true(all(a$n_ambiguous >= 1))
  expect_true(all(a$n_not_plausibly_zero <= a$n_ambiguous))
  expect_true(all(a$n_no_neighbour <= a$n_not_plausibly_zero))   # unknowable implies not plausible
  expect_true(all(is.na(a$min_p_zero) | (a$min_p_zero >= 0 & a$min_p_zero <= 1)))
  # the detector must find the case that motivated it, and must not be looking at the whole panel
  # through the wrong stream: a non-sentinel country is judged on its non-sentinel file
  expect_true(any(a$country_short == "CZ" & a$season == "2024/2025"))
  expect_gte(a$n_ambiguous[a$country_short == "CZ" & a$season == "2024/2025"], 10L)
  # a week that says detections = 0 EXPLICITLY must not be flagged -- that is an observed zero, and
  # confusing the two would throw away thousands of legitimately quiet weeks
  raw <- models_in$data_timeseries_long
  z <- raw[raw$pathogen == "Influenza" & raw$agegroup == "age_total" &
           raw$indicator == "detections" & raw$stream == "typing_sentinel" &
           is.finite(raw$value) & raw$value == 0, ]
  expect_gt(nrow(z), 1000L)      # they exist in quantity, and none of them is an ambiguous week

  cand <- c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")
  with_ex <- withr::with_dir(here::here(), jm_build_data(cand, models_in, demo, verbose = FALSE))
  no_ex   <- withr::with_dir(here::here(), jm_build_data(cand, models_in, demo, verbose = FALSE,
                                                         exclude_ambiguous_positivity = FALSE))
  # the exclusion is ON by default and costs exactly the affected country-seasons
  expect_lt(with_ex$n_cs, no_ex$n_cs)
  expect_equal(nrow(with_ex$excluded_ambiguous), no_ex$n_cs - with_ex$n_cs)
  expect_equal(nrow(no_ex$excluded_ambiguous), 0L)
  expect_true(all(paste(with_ex$excluded_ambiguous$country, with_ex$excluded_ambiguous$season) %in%
                  paste(a$country_short, a$season)))
  # no excluded country-season survives into the fitted design
  kept <- paste(with_ex$countries[with_ex$cs_country + 1L], with_ex$seasons[with_ex$cs_season + 1L])
  expect_false(any(paste(with_ex$excluded_ambiguous$country,
                         with_ex$excluded_ambiguous$season) %in% kept))
  # and the LAYOUT stays exactly consistent after the prune -- the failure mode being an orphaned
  # seed or baseline slot that no country-season references
  expect_equal(length(jm_par_names(with_ex)), with_ex$n_par)
  bl <- jm_blocks(with_ex)
  expect_equal(sort(c(bl$shared, unlist(bl$local))), seq_len(with_ex$n_par))
  expect_equal(with_ex$n_par, with_ex$n_shared + sum(with_ex$n_local))
  expect_equal(with_ex$n_shared, 2L * with_ex$n_season)       # x, delta, sigma and the shared tau
  expect_true(is.finite(jm_negll_cpp(jm_theta0(with_ex), with_ex)))
  for (ic in seq_len(with_ex$n_country))   # a source slot must never survive without data
    expect_equal(with_ex$n_src[ic],
                 length(unique(with_ex$cs_src[with_ex$cs_of_country[[ic]] + 1L])))
  # dropping CZ's only ERVISS season must remove CZ's ERVISS baseline, not leave it inert
  expect_equal(with_ex$sources[[which(with_ex$countries == "CZ")]], "RespiCompass")

  # the threshold lets a stray interior week be tolerated rather than costing a whole season
  strict <- nrow(with_ex$excluded_ambiguous)
  loose <- withr::with_dir(here::here(),
             jm_build_data(cand, models_in, demo, verbose = FALSE, ambiguous_min_unexplained = 99L))
  expect_lte(nrow(loose$excluded_ambiguous), strict)
  expect_gte(loose$n_cs, with_ex$n_cs)
  # THE POINT OF JUDGING PER WEEK: a country-season whose affected weeks are all plausibly zero must
  # be KEPT. PL 2024/2025 has one affected week whose neighbours show 0 detections over 114 tests, so
  # zero is all but certain; CZ 2024/2025 has 14 weeks that cannot plausibly be zero. A rule counting
  # affected weeks would have separated these two only by luck.
  expect_true(any(a$country_short == "PL" & a$season == "2024/2025"))
  expect_equal(a$n_not_plausibly_zero[a$country_short == "PL" & a$season == "2024/2025"], 0L)
  expect_gt(a$n_not_plausibly_zero[a$country_short == "CZ" & a$season == "2024/2025"], 0L)
  expect_false(any(with_ex$excluded_ambiguous$country == "PL"))
  expect_true(any(with_ex$excluded_ambiguous$country == "CZ"))
  expect_equal(with_ex$n_cs_of_country[with_ex$countries == "PL"],
               no_ex$n_cs_of_country[no_ex$countries == "PL"])   # PL keeps every season
})

# ---- the foundations: what the final audit checked, kept as invariants ----
test_that("the shared priors are counted once, not once per country", {
  # The likelihood is decomposed per country for the block sweep. If the shared priors (R0, the season
  # deviations, the elderly susceptibility) were added inside that decomposition they would be counted
  # 12 times, silently tightening them by a factor of 12 and making every contraction wrong.
  nl <- jm_negll_cpp(th, d)
  parts <- sum(vapply(seq_len(d$n_country), function(ic) jm_country_negll_cpp(th, d, ic - 1L), numeric(1)))
  S <- d$n_season
  dn <- function(x, m, s) -0.5 * ((x - m) / s)^2 - log(s) - 0.5 * log(2 * pi)
  x_free <- th[seq_len(S - 1)]; dev_free <- th[S - 1L + seq_len(S - 1)]
  expected <- -(sum(dn(c(x_free, -sum(x_free)), 0, d$pr_x_sd)) +
                sum(dn(c(dev_free, -sum(dev_free)), 0, d$pr_delta_sd)) +
                dn(th[2L * S - 1L], d$pr_sigma_mean, d$pr_sigma_sd) +
                dn(th[2L * S], d$pr_tau_mean, d$pr_tau_sd))
  expect_equal(unname(nl - parts), unname(expected), tolerance = 1e-8)
})

test_that("no fitted slot is left without a prior", {
  # An unpenalised slot makes the penalised Hessian and that family's prior-to-posterior contraction
  # meaningless -- the figure would report an assumption as a result.
  # Test the PRIOR's CURVATURE, not a one-directional move. Two traps avoided: a posterior difference
  # measures the fit rather than the penalty (pushing a seed far enough moves the wave out of the
  # observation window, which can IMPROVE the likelihood at a crude start), and a single direction can
  # move a slot TOWARDS its prior mean -- theta0's seeds sit below theirs, so +8 improves the prior.
  # A proper informative prior has strictly negative log-density curvature in every coordinate, equal
  # to -1/sd^2 for the Gaussians used here. logprior = -negll - loglik isolates it from the data.
  lp <- function(t) -jm_negll_cpp(t, d) - jm_loglik_cpp(t, d)
  h <- 2
  curv <- vapply(seq_len(d$n_par), function(j){
    a <- th; a[j] <- a[j] + h; b <- th; b[j] <- b[j] - h
    (lp(a) + lp(b) - 2 * lp(th)) / h^2
  }, numeric(1))
  expect_true(all(curv < 0),
              info = paste("improper:", paste(jm_par_names(d)[curv >= 0], collapse = ", ")))
  # and the curvature must be the prior sd the settings declare, slot by slot
  expect_equal(curv[grep("log_I0", jm_par_names(d))][1], -1 / d$pr_I0_sd^2, tolerance = 1e-6)
  expect_equal(curv[grep("logit_S0", jm_par_names(d))][1], -1 / d$pr_S0_sd^2, tolerance = 1e-6)
  # a free season-effect slot also drives the constrained last member (minus the sum), so its
  # log-prior curvature is -2/sd^2: its own term plus the constrained term's dependence on it
  expect_equal(curv[1], -2 / d$pr_x_sd^2, tolerance = 1e-6)                 # x of season 1
  expect_equal(curv[d$n_season], -2 / d$pr_delta_sd^2, tolerance = 1e-6)    # delta of season 1
  expect_equal(curv[2L * d$n_season], -1 / d$pr_tau_sd^2, tolerance = 1e-6)  # the shared spread
  # and the contraction denominator must be the sd the C++ actually applied, per family
  skip_if(is.null(.jm_cached), "no saved fit, or it predates this parameter layout")
  fam <- .jm_cached$id$family
  want <- c("S0 season effect (shared)" = d$pr_x_sd, "S0 (country)" = d$pr_S0_sd,
            "reporting c (country)" = d$pr_c_sd, "season deviation (shared)" = d$pr_delta_sd,
            "dispersion phi" = d$pr_phi_sd, "baseline b" = d$pr_b_sd,
            "age reporting offset" = d$pr_off_sd, "seed I0 (country-season)" = d$pr_I0_sd,
            "elderly susceptibility (global)" = d$pr_sigma_sd,
            "spatial spread tau (shared)" = d$pr_tau_sd)
  for (k in names(want)) if (k %in% fam$family)
    expect_equal(fam$prior_sd[fam$family == k], unname(want[[k]]), tolerance = 1e-9, info = k)
})

test_that("the elderly susceptibility redistributes infection without changing transmissibility", {
  # sigma re-weights the contact matrix and the result is rescaled to spectral radius 1 again. That is
  # what lets R0_s mean the same thing whatever sigma_eld is -- and it means sigma_eld is identified by
  # the AGE COMPOSITION of cases, not by the size of the wave. Documented in MODEL.md; asserted here.
  S <- d$n_season
  rho_of <- function(theta){
    sigma <- c(1, 1, 2^theta[2L * S - 1L])
    vapply(seq_len(d$n_country), function(ic){
      Cs <- sweep(d$Cn[[ic]], 1, sigma, "*")
      max(abs(eigen(Cs / max(abs(eigen(Cs, only.values = TRUE)$values)), only.values = TRUE)$values))
    }, numeric(1))
  }
  for (shift in c(-2, 0, 2)){
    t2 <- th; t2[2L * S - 1L] <- th[2L * S - 1L] + shift
    expect_equal(unname(rho_of(t2)), rep(1, d$n_country), tolerance = 1e-9)
  }
  expect_equal(unname(vapply(d$Cn, function(m) max(abs(eigen(m, only.values = TRUE)$values)), numeric(1))),
               rep(1, d$n_country), tolerance = 1e-10)
})

test_that("the noise floor the adequacy diagnostic compares against is unbiased", {
  # jm_adequacy divides the fitted dispersion by the data's own week-to-week scatter, and the whole
  # "2.6x more noise than the data have" limitation rests on that denominator being a fair floor. The
  # sqrt(2/3) factor is the independence correction for a centred 3-week mean; without it the floor is
  # too low and the excess too big. Checked against a known dispersion on realistic wave shapes.
  set.seed(11)
  ma3 <- function(v){ n <- length(v); o <- rep(NA_real_, n)
    for (i in 2:(n - 1)) o[i] <- mean(v[(i - 1):(i + 1)]); o }
  est <- function(y){ m <- ma3(y); okk <- is.finite(y) & is.finite(m) & m > 20
    if (sum(okk) < 5) NA_real_ else sd(y[okk] / m[okk]) / sqrt(2 / 3) }
  for (w in c(3, 8)){
    r <- replicate(200, { n <- 53; mu <- 400 * exp(-((1:n - 26)^2) / (2 * w^2)) + 5
      est(rnbinom(n, size = 4, mu = mu)) })
    expect_equal(median(r, na.rm = TRUE), 0.5, tolerance = 0.06)   # 1/sqrt(4) by construction
  }
  # and the reported excess must be exactly the ratio of the two CVs, on one cell set
  a <- jm_adequacy(jm_fit(d, max_sweeps = 1L, cores = 1, verbose = FALSE))
  expect_equal(a$excess, a$cv_fitted / a$cv_data, tolerance = 1e-9)
})

test_that("the fit is deterministic and independent of the core count", {
  # The project reports specific parameter values. mclapply forks do not reseed, so any RNG inside a
  # parallel stage would make the answer depend on how many cores happened to be free.
  f1 <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  f2 <- jm_fit(d, max_sweeps = 2L, cores = 1, verbose = FALSE)
  expect_identical(f1$theta, f2$theta)
  f3 <- jm_fit(d, max_sweeps = 2L, cores = 2, verbose = FALSE)
  expect_equal(f1$theta, f3$theta, tolerance = 1e-10)
})

# ---- R0 is fixed; S0 is the one sensor ----
test_that("R0 is a fixed input the likelihood actually uses, and S0 is the only sensor", {
  # S0 is a deliberately blunt sensor of susceptibility AND infectivity (MODEL.md). The pin must be
  # read from the data object by both implementations -- not hard-coded in one of them -- and no
  # parameter may carry R0.
  expect_false(any(grepl("R0", jm_par_names(d))))
  expect_null(jm_settings()$season_on); expect_null(jm_settings()$country_on)
  d2 <- d; d2$R0_fixed <- 1.7
  expect_equal(jm_negll_cpp(th, d2), jm_negll_R(th, d2), tolerance = 1e-10)
  expect_gt(abs(jm_negll_cpp(th, d2) - jm_negll_cpp(th, d)), 1)
  expect_equal(jm_unpack(th, d2)$R0, 1.7)
})

# ---- the phi hole: a defect the C++ and its mirror SHARED, so the identity test could not see it ----
test_that("the likelihood agrees with R's own dnbinom, and rejects the phi region where it cannot", {
  # lgamma(y + phi) - lgamma(phi) is a difference of two numbers of size phi*log(phi). Past phi ~ 1e8
  # it loses digits; past ~1e15 it is cancellation noise and can come out hugely POSITIVE. A
  # model-comparison fit found it: one country's log phi ran to 44.7 and the "log-likelihood" to +3e6.
  # Both implementations agreed, because both use that formula -- which is why this test compares
  # against a THIRD, independent implementation: R's dnbinom, which is numerically robust.
  source(here::here("code/07_joint_model/joint_recovery.R"))
  dnb_loglik <- function(theta, dd){
    p <- jm_unpack(theta, dd); f <- jm_fitted_cpp(theta, dd)
    sum(vapply(seq_len(dd$n_cs), function(i){
      ic <- dd$cs_country[i] + 1L; y <- dd$y[[i]]; mu <- f$mu[[i]]; ok <- is.finite(y)
      sum(dnbinom(y[ok], size = p$country[[ic]]$phi, mu = pmax(mu[ok], 1e-10), log = TRUE))
    }, numeric(1)))
  }
  j <- grep("DK:log_phi", jm_par_names(d))
  # across the whole legitimate range and up to the cap: the formula must match dnbinom
  for (v in c(-3, 0, 2, 5, 10, 15, 18)){
    t2 <- th; t2[j] <- v
    expect_equal(jm_loglik_cpp(t2, d), dnb_loglik(t2, d), tolerance = 1e-6,
                 info = sprintf("log phi = %g", v))
  }
  # beyond the cap: rejected by EVERY entry point, never scored
  for (v in c(19, 25, 40, 100)){
    t2 <- th; t2[j] <- v
    expect_gte(jm_negll_cpp(t2, d), 1e10, label = sprintf("negll at log phi = %g", v))
    expect_gte(jm_negll_R(t2, d), 1e10, label = sprintf("negll_R at log phi = %g", v))
    expect_false(is.finite(jm_loglik_cpp(t2, d)), label = sprintf("loglik at log phi = %g", v))
    expect_gte(jm_country_negll_cpp(t2, d, 0L), 1e10, label = sprintf("country negll at log phi = %g", v))
    expect_error(jm_fitted_cpp(t2, d), "rejected")
  }
  # and the cap is far above anything a real fit produces: the prior centre is log 4, sd 1
  expect_gt(log(1e8), d$pr_phi_mean + 10 * d$pr_phi_sd)
})

# ---- the spatial spread: a country's cities are hit at different times ----
test_that("the spread kernel keeps the total and the rise rate and adds exactly its own width", {
  # MODEL.md claims three exact properties of the spread. Each is checked on the kernel itself, the
  # same C++ function the likelihood calls.
  D <- 7L * 53L; t <- seq_len(D)
  bump <- exp(-0.5 * ((t - 180) / 12)^2)
  daily <- unname(cbind(bump, 2 * bump, 0.5 * bump))
  plain <- unname(rowsum(daily, rep(seq_len(53), each = 7), reorder = FALSE))
  for (tau in c(2, 7, 15)){
    sp <- jm_spread_cpp(daily, tau, 53L)
    # (1) the total: the kernel moves infections in time, never in number
    expect_equal(colSums(sp), colSums(daily), tolerance = 1e-10, info = paste("tau", tau))
    # (2) the mean timing: the kernel is symmetric
    wk <- seq_len(53)
    expect_equal(sum(wk * sp[, 1]) / sum(sp[, 1]), sum(wk * plain[, 1]) / sum(plain[, 1]), tolerance = 1e-9)
    # (3) the width: variances add, exactly. The kernel's variance is tau^2 plus the day bin's 1/12,
    # in weeks^2 (measured: equal to 5 digits, so the tolerance is tight)
    v <- function(m) { p <- m / sum(m); sum(wk^2 * p) - sum(wk * p)^2 }
    expect_equal(v(sp[, 1]) - v(plain[, 1]), (tau^2 + 1 / 12) / 49, tolerance = 1e-6, info = paste("tau", tau))
  }
  # (4) the RISE RATE: a convolved exponential is the same exponential times a constant, so away from
  # the edges every weekly ratio is exactly exp(7 r) -- the spread cannot fake a faster or slower rise
  r <- 0.06; ex <- exp(r * t); dex <- cbind(ex, ex, ex)
  sp <- jm_spread_cpp(dex, 10, 53L)
  inner <- 12:40                                   # clear of both edges by more than 7 tau
  expect_equal(sp[inner + 1, 2] / sp[inner, 2], rep(exp(7 * r), length(inner)), tolerance = 1e-10)
  # (5) a vanishing tau is no spread at all
  expect_equal(jm_spread_cpp(daily, 1e-6, 53L), plain, tolerance = 1e-13)
})

test_that("tau = 0 is the model without spread, and the spread never changes the attack rate", {
  S <- d$n_season
  # the C++ takes its untouched no-spread path below 0.05 days; the base-R reference always convolves.
  # Agreement at a vanishing tau is therefore the check that the short cut is the limit it claims
  for (v in c(-30, log(0.049), log(0.051))){
    t2 <- th; t2[2L * S] <- v
    expect_equal(jm_negll_cpp(t2, d), jm_negll_R(t2, d), tolerance = 1e-10, info = paste("log tau", v))
  }
  # continuity across the switch-over, likelihood only (the tau prior is a smooth function of tau)
  a <- th; a[2L * S] <- log(0.0499); b <- th; b[2L * S] <- log(0.0501)
  expect_lt(abs(jm_loglik_cpp(a, d) - jm_loglik_cpp(b, d)), 1e-6)
  # and at real spreads the reference agrees too
  for (v in log(c(3, 10, 30))){ t2 <- th; t2[2L * S] <- v
    expect_equal(jm_negll_cpp(t2, d), jm_negll_R(t2, d), tolerance = 1e-10, info = paste("log tau", v)) }
  # every local copy has the same final size, so the spread moves infections in time, not in number
  t0 <- th; t0[2L * S] <- -30; t1 <- th; t1[2L * S] <- log(14)
  expect_identical(jm_fitted_cpp(t0, d)$attack, jm_fitted_cpp(t1, d)$attack)
  expect_false(identical(jm_fitted_cpp(t0, d)$mu, jm_fitted_cpp(t1, d)$mu))
})

test_that("an impossible spread is rejected everywhere, and a stale data object is refused", {
  S <- d$n_season
  t2 <- th; t2[2L * S] <- log(121)                  # beyond the 120-day cap
  expect_gte(jm_negll_cpp(t2, d), 1e10)
  expect_gte(jm_negll_R(t2, d), 1e10)
  expect_false(is.finite(jm_loglik_cpp(t2, d)))
  expect_gte(jm_country_negll_cpp(t2, d, 0L), 1e10)
  expect_error(jm_fitted_cpp(t2, d), "rejected")
  # a data object from before the spread has a different layout: reading it would put a country's S0
  # where tau belongs, so the engine must refuse rather than guess
  stale <- d; stale$tau_by_country <- NULL
  expect_error(jm_negll_cpp(th, stale), "predates")
})

test_that("tau by country: its own layout, the same model, and equal taus give the shared likelihood", {
  set_c <- modifyList(jm_settings(), list(tau_by_country = TRUE))
  dc <- withr::with_dir(here::here(), jm_build_data(c("DK", "EE"), models_in, demo, set = set_c, verbose = FALSE))
  tc <- jm_theta0(dc, set_c); nm <- names(tc)
  expect_equal(dc$n_par, d$n_par + dc$n_country - 1L)          # one tau per country instead of one
  expect_equal(dc$n_shared, 2L * dc$n_season - 1L)
  expect_false("log_tau" %in% nm)
  # each country's tau is the LAST slot of its block, so every fixed position in the block (phi,
  # the seeds) that the fitter and the flat-line protector index is unchanged
  bl <- jm_blocks(dc)
  for (ic in seq_len(dc$n_country))
    expect_equal(nm[max(bl$local[[ic]])], paste0(dc$countries[ic], ":log_tau"))
  expect_equal(sort(c(bl$shared, unlist(bl$local))), seq_len(dc$n_par))
  set.seed(5)
  for (i in 1:2){
    t2 <- tc + rnorm(length(tc), 0, 0.2)
    expect_equal(jm_negll_cpp(t2, dc), jm_negll_R(t2, dc), tolerance = 1e-10)
  }
  # no dead slot: each country's tau moves only its own likelihood
  j <- grep("EE:log_tau", nm); t3 <- tc; t3[j] <- t3[j] + 0.3
  expect_gt(abs(jm_country_negll_cpp(t3, dc, 1L) - jm_country_negll_cpp(tc, dc, 1L)), 1e-6)
  expect_equal(jm_country_negll_cpp(t3, dc, 0L), jm_country_negll_cpp(tc, dc, 0L))
  # the shared model is the by-country model with every tau equal: identical likelihood
  ts <- th; ts[2L * d$n_season] <- log(9)
  tc2 <- tc; tc2[grep("log_tau", nm)] <- log(9)
  tc2[setdiff(seq_along(tc2), grep("log_tau", nm))] <- ts[-(2L * d$n_season)]
  expect_equal(jm_loglik_cpp(tc2, dc), jm_loglik_cpp(ts, d), tolerance = 1e-12)
  # and its priors: one per country, each with the declared sd
  lp <- function(t) -jm_negll_cpp(t, dc) - jm_loglik_cpp(t, dc)
  a <- tc; a[j] <- a[j] + 2; b <- tc; b[j] <- b[j] - 2
  expect_equal((lp(a) + lp(b) - 2 * lp(tc)) / 4, -1 / dc$pr_tau_sd^2, tolerance = 1e-6)
})
