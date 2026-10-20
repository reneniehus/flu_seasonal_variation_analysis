# Contact-matrix rules of the compartmental model (code/06_comp_model/ASSUMPTIONS.md B3, C3).
# These pin the R0 calibration decision of 2026-09: the force of infection must realise EXACTLY
# the intended R0 at full susceptibility, given the age mixing -- and they document, as a number, the
# artefact of the previous (Stan) convention that the decision replaced.

source(here::here("code/06_comp_model/contact_matrix.R"))
source(here::here("code/06_comp_model/comp_model_settings.R"))

# a small assortative test matrix with unequal populations (rows = contacting group)
C4 <- matrix(c(6, 2, 3, 0.5,
               1, 8, 4, 0.5,
               0.4, 1.2, 9, 1.5,
               0.3, 0.6, 3, 4), 4, 4, byrow = TRUE)
N4 <- c(0.5, 1.2, 6.5, 1.8) * 1e6
C4 <- (C4 * N4 + t(C4 * N4)) / (2 * N4)          # impose reciprocity of directed ends, as the data have

test_that("collapsing groups conserves directed contact-ends and keeps reciprocity", {
  cg <- collapse_contact_groups(C4, N4, list(young = 1:2, medium = 3, elderly = 4))
  expect_equal(sum(cg$C * cg$N), sum(C4 * N4))                    # total ends unchanged
  E3 <- cg$C * cg$N
  expect_equal(E3, t(E3), tolerance = 1e-12)                       # ends a->j == ends j->a
  expect_equal(unname(cg$N), c(N4[1] + N4[2], N4[3], N4[4]))
})

test_that("spectral-radius scaling realises the intended R0 exactly, coordinate-free", {
  Cn <- contact_matrix_normalised(C4, N4, "spectral_radius")
  expect_equal(spectral_radius(Cn), 1, tolerance = 1e-10)
  expect_equal(realised_R0(Cn, 1.5), 1.5, tolerance = 1e-10)
  D <- diag(N4)                                                    # count coordinates: D C D^-1
  expect_equal(spectral_radius(D %*% Cn %*% solve(D)), 1, tolerance = 1e-10)
})

test_that("the Stan convention (beta/cbar * C) overshoots R0 for an assortative matrix", {
  # rho(C) >= cbar with equality only under proportionate mixing -> realised R0 > 1.5
  Cs <- contact_matrix_normalised(C4, N4, "stan_cbar")
  expect_gt(realised_R0(Cs, 1.5), 1.5)
  # proportionate mixing (every group has the same contact vector scaled by population shares) is the
  # boundary case where the two conventions agree
  Cp <- matrix(rep(N4 / sum(N4), each = 4), 4, 4) * 12
  expect_equal(realised_R0(contact_matrix_normalised(Cp, N4, "stan_cbar"), 1.5), 1.5, tolerance = 1e-10)
})

test_that("on the real country matrices the old convention realised R0 = 1.58-1.73 (documented artefact)", {
  skip_if_not(file.exists(here::here("output/models_in.rds")))
  skip_if_not(file.exists(here::here("output/demography_respicast.Rdata")))
  ct <- readRDS(here::here("output/models_in.rds"))$contacts
  load(here::here("output/demography_respicast.Rdata"))           # -> obj$population_pyramid
  g4 <- c("0-4", "5-14", "15-64", "65+")
  settings <- comp_model_settings()
  r0 <- vapply(setdiff(names(ct), "EU"), function(cc){
    pyr <- obj$population_pyramid[obj$population_pyramid$country == cc, ]
    N   <- pyr$population[match(g4, pyr$age_group)]
    cm  <- model_contact_matrix(ct[[cc]], N, settings)              # the model's 3-group scaled matrix
    stan <- contact_matrix_normalised(collapse_contact_groups(ct[[cc]], N, list(1:2, 3, 4))$C,
                                      cm$N, "stan_cbar")
    c(model = realised_R0(cm$C, 1.5), stan = realised_R0(stan, 1.5))
  }, numeric(2))
  expect_true(all(abs(r0["model", ] - 1.5) < 1e-8))                # the fix: exactly 1.5 everywhere
  expect_true(all(r0["stan", ] > 1.55 & r0["stan", ] < 1.80))      # the artefact, as documented
  expect_gt(diff(range(r0["stan", ])), 0.1)                        # ...and it varied by country
})
