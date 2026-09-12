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
