# comp_model_settings.R
#
# EVERY fixed input, prior and structural switch of the compartmental model, in one place, with its
# provenance. The prose version is code/06_comp_model/ASSUMPTIONS.md (item tags A1..I2 refer to it);
# the two files must agree -- change both or neither.
#
# Provenance tags: [owner] decided 2026-09 · [stan] the parked Stan model · [notes] the owner's old
# settings notes · [data] forced by data availability · [proposal] my recommendation.

comp_model_settings = function(){
  p = list()

  # ---- |-age structure (B1-B4) ----
  p$age_groups      = c(young = "0-14", medium = "15-64", elderly = "65+")   # [owner][data]
  p$age_bands_data  = list(young = c("age_00_04", "age_05_14"),               # the ERVISS/RespiCompass bands
                           medium = "age_15_64", elderly = "age_65_99")       # that each model group merges (B2)
  p$contact_bands   = list(young = 1:3, medium = 4:13, elderly = 14:17)       # 17 five-year Prem bands (B3)

  # ---- |-transmission (C1-C6) ----
  p$gamma_per_day   = 0.2777778   # [notes] recovery rate; infectious period 3.6 d (NB: susc_* methods use 3 d)
  p$R0_reference    = 1.5         # [notes] the transmissibility every season is pulled towards (C4)
  p$R0_by_season    = TRUE        # [owner] R0_s free per season, shared across countries (C4)
  p$prior_logR0_sd  = 0.05        # [proposal] log(R0_s) ~ N(log 1.5, 0.05): strong pull, ~10% at 2 sd
  p$ve_spread       = 0.20        # [notes] vaccinated infectiousness x (1 - ve_spread) (C6)
  p$dt_days         = 1           # [stan] daily forward-Euler grid (C5); the Stan runner used n_daily_time_steps = 1
  # C3 -- THE R0 CALIBRATION RULE (read ASSUMPTIONS.md C3 for the numbers):
  #   the contact matrix enters the force of infection ONLY after division by its dominant eigenvalue,
  #   Cn = C / rho(C), so that the next-generation matrix at full susceptibility, (beta/gamma) * Cn, has
  #   spectral radius exactly R0_s. The Stan model's (beta / cbar) * C convention realised R0 = 1.58-1.73
  #   (median 1.65) instead of 1.5, varying by country. contact_matrix_normalised() in the model code
  #   implements this rule and asserts it (rho(Cn) == 1 within 1e-10).
  p$contact_scaling = "spectral_radius"   # [owner][proposal]; "stan_cbar" is kept ONLY to reproduce the old behaviour in tests

  # ---- |-vaccination (D1-D5) ----
  p$ve_inf          = 0.25        # [notes] VE against infection given exposure (susceptibility)
  p$ve_ili_cond_inf = 0.20        # [notes] VE against ILI given infection (severity)
  p$ve_source       = "fixed"     # [owner] "fixed" = the notes values; "csv_by_season" = season VE vs dominant subtype (D3)
  p$vax_pulse_monthday = "-10-01" # [stan][owner] the season's 65+ coverage applied as ONE pulse on 1 October (D2)
  p$vax_groups      = "elderly"   # [stan][data] only the elderly are vaccinated; young/medium coverage = 0 (D2)

  # ---- |-seasons and initial conditions (E1-E4) ----
  # THE TWO-WAY DESIGN (owner, 2026-09-10): S0 is per COUNTRY and shared across seasons; R0 is per
  # SEASON and shared across countries. The rise rate of country-season (c,s) is gamma*(R0_s*S0_c - 1),
  # a season factor times a country factor, so both are identified (per country once S0_c is shared
  # across its seasons; across countries by pooling R0_s).
  p$S0_by_season    = FALSE       # [owner] one S0 per country, shared across seasons (E2)
  p$S0_by_age       = FALSE       # [proposal] ... and shared across age groups (E2)
  p$I0_fraction     = 1e-5        # [proposal] fixed seed in every age group on season day 1 (E3)
  p$season_start_monthday = "-08-01"   # [data] as the panel
  p$reset_each_season = TRUE      # [stan] compartments reset at every season start; no immunity carry-over (E1)

  # ---- |-observation model (F1-F5) ----
  p$detection_age_invariant = TRUE   # [owner] one positivity and one reporting proportion for all ages (F1, F2)
  p$c_by_season     = FALSE       # [open]  reporting proportion c per country, shared across seasons (F2)
  p$rate_per        = 1e5         # [stan]  ILI+ rates are per 100 000 of the age group; counts = rate * pop / 1e5 (F3)
  p$use_cum_burden_term = FALSE   # [owner] no separate cumulative-burden likelihood, no likelihood weights (F4)
  p$obs_weights     = NULL        # [owner] weight_obs_epi dropped (F4)

  # ---- |-process noise / EKF (G1-G3) ----
  p$noise_on        = "log_I"     # [proposal] multiplicative process noise on the infected fractions (G2)
  p$prior_logq      = c(mean = log(0.05), sd = 1)   # [open] regularise q small
  p$p0_frac         = 0.05        # [stan-r] initial-state sd as a fraction of S0 / I0 (G3)

  # ---- |-priors as penalties (H1-H3) ----
  p$prior_logitS0   = c(mean = qlogis(0.75), sd = 1)   # [proposal]
  p$prior_logphi    = c(mean = log(15), sd = 0.8)      # [proposal]

  # ---- |-inference (I1-I2) ----
  # FEASIBILITY NOTE. State per country-season = 3 ages x (S_u,I_u,S_v,I_v) + 3 weekly-incidence
  # accumulators = 15 (R_u, R_v are implied by the within-group sum). A weekly EKF step costs one
  # 7-day integration plus a Jacobian (analytical, or 15 finite-difference integrations). In base R the
  # interpreter overhead dominates: ~30-50 ms per country-season likelihood, ~5-8 s for the 166
  # country-seasons of a joint fit, and a joint MAP over ~230 parameters (166 S0 + 8 R0_s + 25 c + 25 b
  # + phi, q) with numerical gradients needs ~460 likelihood calls per gradient -- ~1 h per BFGS step in
  # base R, i.e. NOT feasible. In Rcpp the same likelihood is ~0.5 ms per country-season, ~0.1 s joint,
  # ~1 min per gradient, hours per joint fit -- feasible, and exact gradients (TMB/autodiff) would make
  # it minutes. So: base R for the reference model and the per-country stage; Rcpp for the joint fit.
  p$stage           = "per_country"   # "per_country" (S0_c + regularised R0_{c,s} per country) | "joint" (R0_s shared across countries)
  p$n_starts        = 4
  p$optim_maxit     = 300
  p$prior_logc      = NULL        # reporting proportion c: unpenalised (data-scale)
  p$prior_logb      = NULL        # baseline b (rate per 100 000): unpenalised

  p
}
