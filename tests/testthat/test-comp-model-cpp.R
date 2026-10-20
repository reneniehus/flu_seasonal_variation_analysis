# The C++ engine must be IDENTICAL to the base-R reference (ASSUMPTIONS.md I2). This test evaluates
# both on real Danish data (when the cached model inputs exist; a synthetic country otherwise) at
# several parameter sets -- including ones with missing weeks, a vaccination pulse, and large process
# noise -- and requires the log-likelihoods and every filtered trajectory to agree to 1e-10 relative.
# If this test fails, one of the two implementations changed without the other: fix that first.

skip_if_not_installed("Rcpp")
source(here::here("code/06_comp_model/contact_matrix.R"))
source(here::here("code/06_comp_model/comp_model_settings.R"))
source(here::here("code/06_comp_model/comp_model_core.R"))
source(here::here("code/06_comp_model/comp_model_cpp.R"))
source(here::here("code/06_comp_model/comp_model_fit.R"))       # cm_b_season, cm_fixed_sigma (the per-season pieces)
withr::with_dir(here::here(), cm_load_cpp())

settings <- comp_model_settings()
C3 <- matrix(c(9, 3, 0.8,  1.5, 7, 1.2,  0.6, 2.5, 3.5), 3, 3, byrow = TRUE); N3 <- c(1.5e6, 6e6, 1.8e6)
C3 <- (C3 * N3 + t(C3 * N3)) / (2 * N3)
f <- cm_fixed(contact_matrix_normalised(C3, N3, "spectral_radius"), N3, settings)

rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8), na.rm = TRUE)

test_that("R and C++ engines agree to 1e-10 on synthetic data (pulse, missing weeks, large q)", {
  set.seed(3)
  mu <- cm_mu(cm_simulate_season(f, 0.72, 1.55, 36, 62, c(0, 0, 0.5))$inc, f, 0.09, 2.5)
  y  <- mu + rnorm(length(mu), 0, sqrt(mu + mu^2 / 20)); y[y < 0] <- 0
  y[c(3, 17, 30), 2] <- NA; y[5, ] <- NA                                     # missing cells and a missing week
  for (par in list(list(S0 = 0.72, R0 = 1.55, c = 0.09, b = 2.5, phi = 20, q = 0.05),
                   list(S0 = 0.85, R0 = 1.40, c = 0.20, b = 0.5, phi = 5,  q = 0.6),
                   list(S0 = 0.60, R0 = 1.70, c = 0.02, b = 8.0, phi = 60, q = 1e-4),
                   list(S0 = 0.75, R0 = 1.50, c = c(0.05, 0.08, 0.16), b = 3, phi = 15, q = 0.1))){   # age-specific reporting
    r <- cm_ekf_season(y, f, par$S0, par$R0, par$c, par$b, par$phi, par$q, 62, c(0, 0, 0.5))
    k <- cm_ekf_season_engine(y, f, par$S0, par$R0, par$c, par$b, par$phi, par$q, 62, c(0, 0, 0.5), engine = "cpp")
    expect_lt(abs(r$loglik - k$loglik) / abs(r$loglik), 1e-10)
    expect_lt(rel(k$mu_pred, r$mu_pred), 1e-10)
    expect_lt(rel(k$I, r$I), 1e-10)
    expect_lt(rel(k$S, r$S), 1e-10)
  }
})

test_that("the C++ deterministic simulator matches the R reference (incidence, end state, attack rates)", {
  for (par in list(list(S0 = 0.8, R0 = 1.5, I0 = 1e-5), list(S0 = 0.6, R0 = 1.7, I0 = 3e-7))){
    r <- cm_simulate_season(f, par$S0, par$R0, 40, 62, c(0, 0, 0.5), I0 = par$I0)
    k <- cm_simulate_season_engine(f, par$S0, par$R0, 40, 62, c(0, 0, 0.5), I0 = par$I0, engine = "cpp")
    expect_lt(rel(k$inc, r$inc), 1e-10)
    expect_lt(rel(as.numeric(k$x_end), r$x_end), 1e-10)
    expect_lt(rel(as.numeric(k$attack), r$attack), 1e-10)
  }
})

test_that("R and C++ engines agree on real Danish data across all seasons", {
  skip_if_not(file.exists(here::here("output/comp_model/fit_DK.rds")))
  fit <- readRDS(here::here("output/comp_model/fit_DK.rds"))
  fd <- cm_fixed(fit$Cn, fit$N, fit$settings); p <- fit$params
  if (!is.null(p$sigma)) fd <- cm_fixed_sigma(fd, p$sigma)            # age-susceptibility profile, if fitted
  for (s in seq_along(fit$seasons)){
    vf <- c(0, 0, fit$vax$coverage[s])
    # per-season pieces of the current parameter layout: c is K x A (season deviation x age offsets),
    # b is per data source, I0 per season -- exactly what cm_negll hands the engines
    cs <- if (is.matrix(p$c)) p$c[s, ] else p$c
    bs <- cm_b_season(p, fit, s)
    I0 <- if (length(p$I0) > 1) p$I0[s] else p$I0
    r <- cm_ekf_season(fit$y[[s]], fd, p$S0, p$R0[s], cs, bs, p$phi, p$q, 62, vf, I0 = I0)
    k <- cm_ekf_season_engine(fit$y[[s]], fd, p$S0, p$R0[s], cs, bs, p$phi, p$q, 62, vf, I0 = I0, engine = "cpp")
    expect_lt(abs(r$loglik - k$loglik) / abs(r$loglik), 1e-10)
    expect_lt(rel(k$mu_pred, r$mu_pred), 1e-10)
  }
})

test_that("the C++ engine is materially faster than the reference (a broken build would not be)", {
  set.seed(4)
  mu <- cm_mu(cm_simulate_season(f, 0.72, 1.55, 40, 62, c(0, 0, 0.5))$inc, f, 0.09, 2.5)
  tR <- system.time(for (i in 1:3) cm_ekf_season(mu, f, 0.72, 1.55, 0.09, 2.5, 20, 0.05, 62, c(0, 0, 0.5)))[["elapsed"]]
  tC <- system.time(for (i in 1:60) cm_ekf_season_engine(mu, f, 0.72, 1.55, 0.09, 2.5, 20, 0.05, 62, c(0, 0, 0.5), engine = "cpp"))[["elapsed"]]
  expect_gt((tR / 3) / (tC / 60), 8)     # measured ~17x on the development machine
})
