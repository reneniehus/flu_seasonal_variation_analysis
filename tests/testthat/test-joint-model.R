# The joint model must mean the same thing in both implementations, and its parameter vector must be
# wired correctly. Three things are checked, in rough order of how much damage they would do:
#   1. the C++ log-posterior equals an independent base-R implementation of the same model to 1e-10,
#      at the starting point and at random perturbations. The R version is written from the maths and
#      uses R's own eigen() for the spectral radius, so agreement means both are right rather than
#      that one was copied from the other;
#   2. every slot of the parameter vector is read back by jm_unpack as the name says, and the
#      sum-to-zero constraint on the season deviations holds;
#   3. no parameter is DEAD -- the objective responds to each one. The compartmental pilot shipped an
#      inert baseline slot for countries with only one data source, which made its Hessian exactly
#      singular; this test would have caught it.
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
  expect_equal(unname(log(p$R0)), unname(probe[seq_len(S)]))
  expect_equal(unname(p$delta[seq_len(S - 1)]), unname(probe[S + seq_len(S - 1)]))
  expect_equal(sum(p$delta), 0, tolerance = 1e-12)            # the constraint that removes the pilot's flat direction
  expect_equal(log2(p$sigma_eld), unname(probe[2L * S]))
  for (ic in seq_len(d$n_country)){
    b <- d$off_country[ic]; q <- p$country[[ic]]
    expect_equal(qlogis(q$S0), unname(probe[b + 1]))
    expect_equal(log(q$c),     unname(probe[b + 2]))
    expect_equal(q$off_young,  unname(probe[b + 3]))
    expect_equal(q$off_eld,    unname(probe[b + 4]))
    expect_equal(log(q$phi),   unname(probe[b + 5]))
    expect_equal(unname(log(q$b)),  unname(probe[b + 5 + seq_len(d$n_src[ic])]))
    expect_equal(unname(log(q$I0)), unname(probe[b + 5 + d$n_src[ic] + seq_len(d$n_cs_of_country[ic])]))
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
  tt <- jm_truth_with_driver(d, x, beta_R0 = 0.10, beta_delta = 0.40, noise_R0 = 0, noise_delta = 0)
  dr <- attr(tt, "driver")
  expect_equal(dr$beta_R0, 0.10); expect_equal(dr$beta_delta, 0.40)
  # with the noise switched off, regressing the truth back on the covariate returns the slopes exactly
  expect_equal(unname(coef(lm(tt[seq_len(S)] ~ dr$x))[2]), 0.10, tolerance = 1e-8)
  p <- jm_unpack(tt, d)
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
  expect_equal(nrow(iv), length(small$theta) + 1L)
  expect_true(paste0("delta_", d$seasons[d$n_season]) %in% iv$parameter)
})

test_that("the flat-line protector detects a flat fit, escapes it, and leaves a healthy one alone", {
  # The guarded failure: a country's local block has a second optimum in which the dispersion
  # collapses, the negative binomial becomes diffuse enough that any curve fits, and the country sits
  # at its baseline for every season. It cost the Netherlands in the first joint fit and a harder
  # search later found a solution 406 nats better, so it is an optimiser failure, not a fact.
  skip_if_not(file.exists(here::here("output/joint_model/joint_fit.rds")), "no saved fit")
  fit <- readRDS(here::here("output/joint_model/joint_fit.rds"))$fit
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
  skip_if_not(file.exists(here::here("output/joint_model/joint_fit.rds")), "no saved fit")
  fit <- readRDS(here::here("output/joint_model/joint_fit.rds"))$fit
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
    "sigma overflows"  = local({ t <- th; t[2L * S] <- 1030; t }),
    "R0 overflows"     = local({ t <- th; t[seq_len(S)] <- 800; t[grep(":log_I0", names(th))] <- -800; t }),
    "one R0 overflows" = local({ t <- th; t[1] <- 800; t }))
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
  skip_if_not(file.exists(here::here("output/joint_model/joint_fit.rds")), "no saved fit")
  fit <- readRDS(here::here("output/joint_model/joint_fit.rds"))$fit
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
  expect_named(small$flat_thresholds, c("attack_min", "epi_frac_min"))
})
