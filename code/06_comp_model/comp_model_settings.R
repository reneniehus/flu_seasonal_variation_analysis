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
  # AGE-SPECIFIC SUSCEPTIBILITY (the alternative to age-specific reporting): relative susceptibility
  # sigma_a for young and elderly (medium = 1), shared across seasons within a country. It enters the
  # force of infection as row scaling of the contact matrix, diag(sigma) %*% Cn, RENORMALISED to
  # spectral radius 1 so R0_s keeps its meaning: sigma redistributes WHO gets infected (attack rates by
  # age), it does not change the overall transmissibility. Discriminates from reporting offsets (F2)
  # through the dynamics -- and through external attack-rate profiles (PHIRST, South Africa).
  p$susc_by_age     = FALSE       # [open] fit log2 sigma_young, log2 sigma_elderly; prior N(0, 1)
  p$prior_log2susc_sd = 1
  p$susc_fixed      = NULL        # or a FIXED profile c(young, 1, elderly) when susc_by_age is FALSE -- the device of
                                  # run_susc_grid.R: one profile SHARED across countries (biology) with reporting
                                  # offsets free per country (surveillance); NULL = no age profile
  p$I0_fraction     = 10^-6.5     # [owner, 2026-09] centre of the per-season seed prior (E3); recentred from 1e-5 after 29% of the
                                  # fitted seeds (5 countries, 35 seasons; median 10^-6.4, range 10^-8.8..10^-3.6) fell below
                                  # the old prior's lower 2.5% bound -- late waves need small seeds
  # SEASON ARRIVAL TIME (E3, evidence: the Danish fit with a fixed seed peaks ~12 weeks too early in
  # every season and the filter has to abandon the SIR). With S0 shared across seasons the seed size is
  # the only smooth handle on WHEN a season arrives: I0_s = I0 * exp(delta_s), delta_s ~ N(0, sd) per
  # season. It shifts the epidemic by delta_s / r days at growth rate r and leaves r = gamma*(R0_s*S0_c-1)
  # -- hence the identification of R0_s and S0_c -- untouched. (The Stan model's disabled i_season term.)
  p$I0_by_season    = TRUE        # [proposal] per-season seed size (arrival time); FALSE = the fixed seed of E3 as first agreed
  p$prior_logI0_sd  = 3           # [owner, 2026-09] wide: 95% band 10^-9.1 .. 10^-3.9 covers every fitted seed so far
  p$season_start_monthday = "-08-01"   # [data] as the panel
  p$reset_each_season = TRUE      # [stan] compartments reset at every season start; no immunity carry-over (E1)

  # ---- |-observation model (F1-F5) ----
  p$detection_age_invariant = TRUE   # [owner] one positivity and one reporting proportion for all ages (F1, F2)
  # Age-specific reporting offsets (the Stan model's prop_ili_age: c_a = c * 2^(off_a), medium = reference,
  # off ~ N(0, 1)). OFF by the owner's decision; the switch exists because the Danish fit under-predicts
  # the elderly ~2x in every season, so the two variants must be comparable side by side.
  p$c_by_age        = TRUE        # [proposal, pending owner] evidence: DK stage-1 fit 118 nats better; medium adults report ~0.5x the young, ~0.35x the elderly
  p$prior_logc_age_sd = 1
  # Per-season reporting deviation (the Stan model's disabled prop_ili_season): c_{c,s} = c_c * exp(delta_s),
  # delta_s ~ N(0, sd). Evidence for needing it: with S0_c shared and R0_s the only season factor, the
  # Danish fit cannot reproduce the 2-4x larger 2015/16 and 2017/18 peaks even under a wide R0 prior
  # (a larger R0 also makes the wave sharper and earlier), so the fit treats them as noise (phi ~ 1).
  # ILI per infection plausibly differs by season/subtype (symptomaticity, care-seeking, positivity).
  p$c_by_season     = TRUE        # [proposal, pending owner] evidence: DK stage-1 fit 89 nats better; deviations 0.5-1.7x track the big/small seasons
  p$prior_logc_season_sd = 0.5    # [proposal] a factor ~2.7 at 2 sd
  p$rate_per        = 1e5         # [stan]  ILI+ rates are per 100 000 of the age group; counts = rate * pop / 1e5 (F3)
  # Off-season baseline b PER DATA SOURCE: RespiCompass ILI+ is exactly 0 in weeks without flu detections,
  # ERVISS-era reconstructions sit at a positive floor; one shared b cannot be both (the Danish fit forced
  # phi ~ 1 to make hundreds of exact zeros plausible under a ~100-count baseline).
  p$b_by_source     = TRUE        # [proposal] log_b per source (RespiCompass, ERVISS) instead of one b
  p$use_cum_burden_term = FALSE   # [owner] no separate cumulative-burden likelihood, no likelihood weights (F4)
  p$obs_weights     = NULL        # [owner] weight_obs_epi dropped (F4)

  # ---- |-process noise / EKF (G1-G3) ----
  p$noise_on        = "log_I"     # [proposal] multiplicative process noise on the infected fractions (G2)
  # q is the weekly multiplicative sd of the infected fractions: 0.05 = 5% wiggle, 0.6 = the filter can
  # nearly reset the state every week. The first Danish fits escaped to q = 0.6 / phi = 1 under a weak
  # prior (log-sd 1), at which point the deterministic SIR no longer mattered to the likelihood -- the
  # 'loose filter masks the model' failure recorded in decisions.md. Hence a firm prior: 95% in [0.04, 0.27].
  p$prior_logq      = c(mean = log(0.1), sd = 0.5)   # [proposal] firm (used only when q_fixed is NULL)
  # q FIXED (owner intent: 'let the trajectory wiggle', not 'let the filter replace the model'). Estimating
  # q from the EKF likelihood let the filter take over (q -> 0.23-0.6) and the backbone parameters drifted
  # away from the stage-1 optimum to an early, oversized deterministic wave. With q fixed at 10% weekly
  # the filter adds modest corrections and the structural parameters stay identified by the data.
  p$q_fixed         = 0.05        # [proposal] 5% weekly; NULL = estimate q with the prior above (see ASSUMPTIONS.md G2)
  # P0 = 0 (G3, decisive Danish evidence): with an initial covariance proportional to the fitted seed,
  # the EKF objective PREFERRED a degenerate mode (negll 6564 vs 9100 at the deterministic optimum):
  # inflating the seed inflated P0, and the filter then carried the wave by state corrections while c, b,
  # phi collapsed. With P0 = 0 the ordering reverses (5573 vs 7366) and the EKF optimum agrees with the
  # deterministic stage in every parameter. The process noise alone provides the wiggle.
  p$p0_frac         = 0           # [proposal] no initial-state uncertainty (the single-population EKF used 0.05)

  # ---- |-priors as penalties (H1-H3) ----
  p$prior_logitS0   = c(mean = qlogis(0.75), sd = 1)   # [proposal]
  p$prior_logphi    = c(mean = log(4), sd = 0.5)       # [proposal] centred on the MEASURED noise: the Danish age-specific series scatter
                                                       # 23% (median) around a 3-week moving average, i.e. phi ~ 3 even for a perfect mean

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
  p$engine          = "cpp"           # "cpp" = the Rcpp port (verified identical to the R reference by test-comp-model-cpp.R; ~17x faster) | "R" = the reference
  p$n_starts        = 4
  p$optim_maxit     = 300
  p$prior_logc      = c(mean = log(0.05), sd = 2)   # reporting proportion: WEAK (a factor ~50 at 2 sd); it exists only to
                                                    # close the 'no epidemic' trap (c -> 1e-200, everything baseline + noise) seen on DK
  p$prior_logb      = NULL        # baseline b (rate per 100 000): unpenalised
  # STAGED FITTING (I1). Stage 1 fits the DETERMINISTIC model (no filter: y ~ N(mu_det, mu + mu^2/phi)),
  # a well-behaved objective that forces the SIR to match timing and shape; stage 2 starts the EKF at
  # that optimum with q regularised. Fitting the EKF from scratch let the filter dominate the likelihood
  # (q -> 0.6, phi -> 1) and the backbone parameters wandered.
  p$two_stage       = TRUE

  p
}
