# Core correctness of the compartmental model engine (code/06_comp_model/comp_model_core.R).
# Four things must hold before any fit is believed: the analytical Jacobian is the derivative of the
# step (the EKF is only as good as this), mass is conserved, the R0 calibration is realised by the
# DYNAMICS (early growth rate = gamma*(R0-1) at full susceptibility), and the EKF collapses onto the
# deterministic model when the process noise and the initial uncertainty vanish.

source(here::here("code/06_comp_model/contact_matrix.R"))
source(here::here("code/06_comp_model/comp_model_settings.R"))
source(here::here("code/06_comp_model/comp_model_core.R"))
source(here::here("code/06_comp_model/comp_model_fit.R"))
source(here::here("code/01_main_supporting/sir_core.R"))          # .fit_multistart

settings <- comp_model_settings()
C3 <- matrix(c(9, 3, 0.8,  1.5, 7, 1.2,  0.6, 2.5, 3.5), 3, 3, byrow = TRUE)
N3 <- c(1.5e6, 6e6, 1.8e6)
C3 <- (C3 * N3 + t(C3 * N3)) / (2 * N3)                          # reciprocity
Cn <- contact_matrix_normalised(C3, N3, "spectral_radius")
f  <- cm_fixed(Cn, N3, settings)

test_that("the analytical day-step Jacobian matches finite differences (with and without the pulse)", {
  # an INTERIOR state (every compartment strictly positive): at a zero compartment the clamp makes the
  # step non-differentiable and a central difference returns half the one-sided derivative
  x <- cm_init(f, 0.7); x[f$L$I_u] <- c(2e-3, 4e-3, 1e-3); x[f$L$S_v] <- c(0.05, 0.1, 0.3); x[f$L$I_v] <- c(1e-4, 2e-4, 5e-4)
  for (vf in list(NULL, c(0, 0, 0.6))){
    st <- cm_day(x, f, beta = 1.5 * f$gamma, vax_frac = vf, jac = TRUE)
    Jfd <- sapply(seq_along(x), function(j){ h <- 1e-7; xp <- x; xp[j] <- xp[j] + h; xm <- x; xm[j] <- xm[j] - h
      (cm_day(xp, f, 1.5 * f$gamma, vf) - cm_day(xm, f, 1.5 * f$gamma, vf)) / (2 * h) })
    expect_equal(st$J, Jfd, tolerance = 1e-6)
  }
})

test_that("the week-map Jacobian matches finite differences", {
  x <- cm_init(f, 0.7); x[f$L$I_u] <- c(1e-3, 2e-3, 5e-4); x[f$L$S_v] <- c(0.02, 0.05, 0.2); x[f$L$I_v] <- c(5e-5, 1e-4, 2e-4)
  st <- cm_week(x, f, 1.6 * f$gamma, day0 = 56, vax_day = 62, vax_frac = c(0, 0, 0.5), jac = TRUE)
  Jfd <- sapply(seq_along(x), function(j){ h <- 1e-6; xp <- x; xp[j] <- xp[j] + h; xm <- x; xm[j] <- xm[j] - h
    (cm_week(xp, f, 1.6 * f$gamma, 56, 62, c(0, 0, 0.5)) - cm_week(xm, f, 1.6 * f$gamma, 56, 62, c(0, 0, 0.5))) / (2 * h) })
  expect_equal(st$J, Jfd, tolerance = 1e-5)
})

test_that("mass is conserved and the vaccination pulse moves exactly the coverage fraction of S_u", {
  sim <- cm_simulate_season(f, S0 = 0.8, R0 = 1.5, n_weeks = 40, vax_day = 62, vax_frac = c(0, 0, 0.6))
  x <- sim$x_end; L <- f$L
  expect_true(all(x[c(L$S_u, L$I_u, L$S_v, L$I_v)] >= 0))
  expect_true(all(x[L$S_u] + x[L$I_u] + x[L$S_v] + x[L$I_v] <= 1 + 1e-12))    # implied R >= 0
  expect_true(all(sim$attack >= 0 & sim$attack <= 0.8))
  # the day of the pulse: S_v jumps by 0.6 * S_u for the elderly only
  x0 <- cm_init(f, 0.8); before <- x0
  for (d in 1:61) before <- cm_day(before, f, 1.5 * f$gamma)
  after <- cm_day(before, f, 1.5 * f$gamma, vax_frac = c(0, 0, 0.6))
  step  <- cm_day(before, f, 1.5 * f$gamma)                                    # same day without the pulse
  expect_equal(after[L$S_v][3], 0.6 * step[L$S_u][3], tolerance = 1e-12)
  expect_equal(after[L$S_v][1:2], c(0, 0))
})

test_that("the R0 calibration is realised by the dynamics: early growth rate = gamma*(R0-1)", {
  # tiny seed, S = 1 - I0 in every group: the total infected grows at gamma*(R0*rho(Cn) - 1) once the
  # dominant eigenvector mode has taken over; with rho(Cn) = 1 that is gamma*(R0-1) exactly
  # (daily Euler carries a small discretisation error of order r^2/2)
  R0 <- 1.5; f2 <- f; f2$I0 <- 1e-8
  x <- cm_init(f2, 1 - 1e-8); tot <- numeric(120)
  for (d in 1:120){ x <- cm_day(x, f2, R0 * f2$gamma); tot[d] <- sum(x[f2$L$I_u] * f2$N) }
  r_hat <- mean(diff(log(tot[60:100])))
  r_true <- log(1 + f2$gamma * (R0 - 1))       # exact per-day multiplier of the Euler linearisation
  expect_equal(r_hat, r_true, tolerance = 5e-3)
  # the Stan convention would have grown faster
  Cs <- contact_matrix_normalised(C3, N3, "stan_cbar"); fs <- cm_fixed(Cs / spectral_radius(Cs), N3, settings)
  expect_gt(spectral_radius(Cs), 1)
})

test_that("with vanishing process noise and initial uncertainty the EKF reproduces the deterministic model", {
  # also with a tiny seed: the filter's post-update clamp must floor at 0 exactly like the deterministic
  # step (a positive floor would re-seed I_v every week and the two would drift apart by floor/I0)
  for (I0 in c(f$I0, 1e-9)){
    sim <- cm_simulate_season(f, 0.75, 1.55, 35, 62, c(0, 0, 0.5), I0 = I0)
    mu  <- cm_mu(sim$inc, f, c = 0.1, b = 2)
    f0 <- f; f0$p0 <- 1e-9
    e <- cm_ekf_season(mu, f0, 0.75, 1.55, 0.1, 2, 20, q = 1e-9, 62, c(0, 0, 0.5), I0 = I0)
    expect_equal(e$mu_pred, mu, tolerance = 1e-9)
    expect_true(is.finite(e$loglik))
  }
})

test_that("parameter recovery on synthetic seasons: S0 and the season R0 ordering come back", {
  # truth chosen so every season completes a wave within 36 weeks at the model's default seed:
  # R_eff = R0*S0 = 1.21 / 1.28 / 1.36 peaks at weeks ~28 / ~21 / ~16 from I0 = 10^-6.5
  set.seed(7)
  S0_true <- 0.85; R0_true <- c(1.42, 1.5, 1.6); c_true <- 0.08; b_true <- 3; phi_true <- 25
  cd <- list(country = "SYN", seasons = c("A", "B", "C"), groups = names(N3), N = N3, Cn = Cn,
             vax_day = 62, vax_frac = list(c(0,0,.5), c(0,0,.5), c(0,0,.5)))
  cd$y <- lapply(seq_along(R0_true), function(s){
    mu <- cm_mu(cm_simulate_season(f, S0_true, R0_true[s], 36, 62, c(0,0,.5))$inc, f, c_true, b_true)
    y  <- mu + rnorm(length(mu), 0, sqrt(mu + mu^2 / phi_true)); y[y < 0] <- 0; y })
  cd$rates <- lapply(cd$y, function(y) sweep(y, 2, N3 / settings$rate_per, "/"))
  cd$n_weeks <- rep(36L, 3)
  cd$vax <- data.frame(season = cd$seasons, coverage = 0.5, provenance = "synthetic")
  cd$source_by_season <- setNames(rep("ERVISS", 3), cd$seasons)      # the per-source baseline needs a source label
  fit <- fit_comp_model(cd, settings, R0_free = TRUE, n_starts = 2, verbose = FALSE)
  expect_equal(fit$convergence, 0L)
  expect_lt(abs(fit$params$S0 - S0_true), 0.06)
  expect_equal(order(fit$params$R0), order(R0_true))
  expect_true(all(abs(fit$params$R0 - R0_true) < 0.05))
})

test_that("start vector, unpack and prior agree on the parameter layout for every switch combination", {
  # a missing slot in any one of the three shifts every later parameter silently (the objective then
  # reads phi off the end of the vector and returns its 1e10 failure value at every start) -- so every
  # combination of the layout switches is exercised: right length, named slots read by cm_unpack,
  # finite objective at the start for both stages, and a prior that responds to exactly its own slots
  set.seed(3); K <- 2L
  cd <- list(country = "SYN", seasons = c("A", "B"), groups = names(N3), N = N3, Cn = Cn, vax_day = 62,
             vax_frac = list(c(0, 0, .5), c(0, 0, .5)), n_weeks = rep(30L, K),
             vax = data.frame(season = c("A", "B"), coverage = 0.5, provenance = "synthetic"),
             source_by_season = c(A = "RespiCompass", B = "ERVISS"))
  cd$y <- lapply(c(1.5, 1.6), function(R0){ mu <- cm_mu(cm_simulate_season(f, 0.85, R0, 30, 62, c(0, 0, .5))$inc, f, 0.08, 3)
    y <- mu + rnorm(length(mu), 0, sqrt(mu + mu^2 / 25)); y[y < 0] <- 0; y })
  cd$rates <- lapply(cd$y, function(y) sweep(y, 2, N3 / settings$rate_per, "/"))
  grid <- expand.grid(c_by_age = c(FALSE, TRUE), c_by_season = c(FALSE, TRUE), susc_by_age = c(FALSE, TRUE),
                      b_by_source = c(FALSE, TRUE), I0_by_season = c(FALSE, TRUE))
  for (i in seq_len(nrow(grid))){
    s <- settings; for (nm in names(grid)) s[[nm]] <- grid[i, nm]
    st <- cm_theta_start(cd, s, R0_free = TRUE); th <- st$base
    n_expected <- 1 + K + (if (s$I0_by_season) K else 0) + 1 + (if (s$c_by_age) 2 else 0) + (if (s$c_by_season) K else 0) +
      (if (s$susc_by_age) 2 else 0) + (if (s$b_by_source) 2 else 1) + 2                 # log_phi, log_q
    expect_equal(length(th), n_expected); expect_equal(length(st$jit), n_expected)
    p <- cm_unpack(th, K, TRUE, s)
    expect_equal(unname(p$phi), unname(exp(th[["log_phi"]])))                      # phi is read from ITS slot
    expect_equal(length(p$b), if (s$b_by_source) 2 else 1)
    expect_true(is.finite(cm_negll(th, cd, f, s, stage = "det")))
    expect_true(is.finite(cm_negll(th, cd, f, s, stage = "ekf")))
    if (s$susc_by_age){ th2 <- th; th2[["log2susc_young"]] <- 1; expect_equal(cm_unpack(th2, K, TRUE, s)$sigma, c(2, 1, 1)) }
    if (s$c_by_age){ th2 <- th; th2[["log2c_elderly"]] <- 1; p2 <- cm_unpack(th2, K, TRUE, s); expect_equal(p2$c_age[3] / p2$c_age[2], 2) }
    # the prior ignores the unpenalised baseline slots and responds to phi
    th3 <- th; th3[grep("^log_b", names(th))] <- th3[grep("^log_b", names(th))] + 1
    expect_equal(cm_logprior(th3, K, TRUE, s), cm_logprior(th, K, TRUE, s))
    th4 <- th; th4[["log_phi"]] <- th4[["log_phi"]] + 1
    expect_false(isTRUE(all.equal(cm_logprior(th4, K, TRUE, s), cm_logprior(th, K, TRUE, s))))
  }
})
