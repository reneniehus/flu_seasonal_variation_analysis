# Compartmental influenza model (`flu_comp_model`) — every assumption, spelled out

Status: agreed with the owner on 2026-09-10 after a full inventory of the parked Stan model
(`stan/SIR_multiseason_age_vax_2.stan`), its recovered input-assembly code (`make_stan_list()` in the
deleted `model_SIR_multiseason.R`, recoverable via `git show 8e1abf1^:code/01_main_supporting/model_SIR_multiseason.R`),
`documentation/documentation.Rmd`, and the data actually in `output/models_in.rds`.
Every item is tagged with its provenance:
`[owner]` = decided by the owner in the 2026-09 discussion · `[stan]` = what the Stan model did ·
`[notes]` = value from the owner's old settings notes · `[data]` = forced by data availability ·
`[proposal]` = my recommendation, adopted unless the owner objects · `[open]` = still to decide.

The machine-readable counterpart is `code/06_comp_model/comp_model_settings.R`; the two must agree.

---

## 0. Purpose

A mechanistic, age- and vaccination-structured SIR that summarises every country-season by an
**initial susceptibility** and lets a **season-level transmissibility** (dominant-subtype driven) be
shared across countries. It is fitted to weekly age-specific ILI+ with an Extended Kalman Filter, so
the latent trajectory may wiggle around the deterministic SIR (process noise). It matches the
assumptions of the previous Stan model except where the owner has explicitly changed them below.

## 1. Population and compartments

- A1 `[stan][owner]` Each age group `a` is modelled as **fractions of that age group's population**;
  compartments S, I, R, each split into unvaccinated (`_u`) and vaccinated (`_v`): six per age group,
  `S_u + I_u + R_u + S_v + I_v + R_v = 1` within each age group.
- A2 `[stan]` Closed population: no births, deaths, migration or ageing within or across seasons.
- A3 `[stan][owner]` No waning of natural or vaccine immunity within a season ("we are not modelling
  waning currently").
- A4 `[stan]` Recovery is a single exponential stage with rate `gamma`; no latent (E) stage.
- A5 `[stan]` Infection confers full protection for the rest of the season (SIR, not SIRS).

## 2. Age structure

- B1 `[owner][data]` Three age groups: **young = 0-14, medium = 15-64, elderly = 65+**. The ERVISS ILI
  bands are exactly 0-4 / 5-14 / 15-64 / 65+ (the RespiCompass ILI+ stream carries the same four), so
  the cut merges the two child bands and is otherwise native to the data.
- B2 `[data]` Merging bands: ILI+ RATES are per 100 000 of the band, so the merged rate is the
  population-weighted mean of the two band rates (equivalently: counts = rate x pop/1e5 add).
- B3 `[data]` Contact matrices come from the Prem et al. synthetic matrices (16 five-year bands,
  extended to 17 with an 80+ band) via `transform_contacts()` in `flu_functions.R`, whose 17->group
  aggregation uses FULL block sums under the directed-contact-ends convention (fixed 2026-08; the
  earlier lower-triangle sum understated within-group mixing ~1.5x). For three groups the same
  aggregation is applied with the cut 0-14 / 15-64 / 65+ (bands 1:3 / 4:13 / 14:17).
- B4 `[data]` Populations per group from the RespiCast pyramid (`data$demography_respicast`).
- B5 `[stan]` Countries without a Prem matrix (Liechtenstein, Norway) use the EU-average matrix.
  (Czechia is now recovered from the workbooks; the Stan runner also used the EU average for CZ.)

## 3. Transmission

- C1 `[stan]` Force of infection on group `a` (per day):
  `lambda_a = beta_s * sum_j Cn[a,j] * ( I_u[j] + (1 - ve_spread) * I_v[j] )`,
  where `I_.[j]` are prevalence FRACTIONS of group `j`, and `Cn` is the contact matrix scaled as in C3.
  Unvaccinated susceptibles are infected at `lambda_a * S_u[a]`, vaccinated at
  `lambda_a * S_v[a] * (1 - ve_inf)`.
- C2 `[stan][notes]` `beta_s = R0_s * gamma`, with `gamma = 0.2777778 / day` (infectious period 3.6 days,
  the owner's old notes; NOTE the single-population susceptibility methods in
  `code/01_main_supporting/methods/` use 3 days = 7/3 per week -- the two are deliberately not
  reconciled here so this model matches the Stan model).
- C3 `[owner][proposal]` **The realised R0 at full susceptibility must equal R0_s exactly, given the age
  mixing.** How: the next-generation matrix at S = 1 (no vaccination) is `K = (beta_s/gamma) * Cn`, so
  R0 = R0_s * rho(Cn) with `rho` the dominant eigenvalue. We therefore use
  `Cn = C / rho(C)` (spectral radius 1). This is coordinate-free (rho is identical for the
  fraction-based and count-based matrices `C` and `D C D^-1`).
  WHY THIS MATTERS -- the Stan model did NOT do this: it used `beta * a_factor[a] * rowNormalised(C)`,
  which algebraically equals `(beta / cbar) * C` with `cbar` the population-weighted mean contacts
  per person. Since rho(C) >= cbar for any assortative matrix, the realised R0 was
  `1.5 * rho(C)/cbar` = **1.58-1.73 across the 28 country matrices (median 1.65), never 1.5**, and it
  differed by country -- a contact-data artefact that would masquerade as between-country
  transmissibility differences. (Computed 2026-09-10 on `output/models_in.rds$contacts`; 3-group
  and 4-group cuts give the same picture.) The `a_factor` device is dropped: `C`'s row sums already
  carry each group's contact activity, and the eigenvalue scaling removes the unreliable absolute
  scale of synthetic matrices (Prem 2017), keeping only the mixing STRUCTURE -- the principled choice
  when R0 is fixed externally (Diekmann, Heesterbeek & Roberts 2010; see decisions.md).
- C4 `[owner]` `R0_s` varies BY SEASON (different dominant subtypes carry different transmission
  potential) and is SHARED ACROSS COUNTRIES; prior: `log(R0_s) ~ N(log 1.5, 0.05)` `[proposal]`,
  i.e. a strong pull to 1.5 with ~10% spread at 2 sd, "to ease fitting few data".
  Identifiability (the two-way design, owner decision 2026-09-10): the initial susceptibility is
  per COUNTRY and shared across seasons (`S0_c`, see E2), so the rise rate of country-season (c, s)
  is `gamma*(R0_s * S0_c - 1)`: a season factor times a country factor. Within one country the
  per-season `R0_{c,s}` IS identified once `S0_c` is shared across its seasons, so the per-country
  stage fits `S0_c` plus regularised per-season `R0_{c,s}` (a diagnostic of how season effects look
  before pooling); the joint stage then pools `R0_s` across countries. (Had `S0` been per
  country-season, `R0_s` and `S0[c,s]` would trade off almost perfectly.)
- C5 `[stan]` Time stepping: forward Euler on a DAILY grid (`dt = 1 day`, `n_sub = 1` as the Stan runner
  used); observations are weekly sums. Sub-daily steps are a settings switch.
- C6 `[stan]` The vaccinated infectious are `(1 - ve_spread)` as infectious as the unvaccinated.
- C7 `[open, 2026-09-10]` AGE-SPECIFIC SUSCEPTIBILITY as the mechanistic alternative to age-specific
  reporting (F2): relative susceptibility `sigma_a` for the young and the elderly (medium = 1; per
  country, shared across seasons; settings `susc_by_age`, prior `log2 sigma ~ N(0, 1)`). It enters the
  force of infection as ROW scaling of the contact matrix, `diag(sigma) %*% Cn`, RENORMALISED to
  spectral radius 1, so `R0_s` keeps its meaning (C3): sigma redistributes WHO gets infected -- the
  attack-rate profile by age -- without changing the overall transmissibility. (The unnormalised
  product would change R0 by rho(diag(sigma) Cn), a pure artefact of the profile.) Contrast with the
  reporting offsets `2^off_a` of F2, which leave the dynamics untouched and rescale the observed
  ILI+ per infection. The two are distinguished by (i) the DYNAMICS -- susceptibility changes each
  group's wave (timing, size, its share of transmission) whereas an offset only rescales it -- and
  (ii) EXTERNAL attack-rate evidence: a reporting-independent prospective cohort (PHIRST, South Africa,
  Cohen et al. 2021, Lancet Glob Health, Figure 2A) finds influenza infection incidence per 100
  person-seasons of 66 (<1 y), 67 (1-4), 51 (5-12), 45 (13-18), 30 (19-44), 30 (45-64), 25 (>= 65):
  young > adult > elderly. The modelled attack-rate profile under each option is compared with this
  ordering in the age experiment (decisions.md, 2026-09). Cross-country logic of that experiment: a
  susceptibility effect should give the SAME direction of over/under-prediction by age in every
  country under age-invariant reporting; reporting differences would vary by surveillance system.
  OUTCOME (12 countries, `run_age_experiment.R` + `run_susc_grid.R`, decisions.md 2026-09): with each
  country's reporting offsets free, one profile shared by all countries is best at young 1.0 and
  elderly 1.7-2.8 (+89 nats total; every country's own best has elderly > 1, and any extra young
  susceptibility loses likelihood in 10 of 12 countries). Reading: the young excess is a reporting
  LEVEL effect; the elderly excess is a shared DYNAMIC effect (elderly ~2x more susceptible per
  contact than the vaccination-adjusted contact matrix implies), which brings the modelled
  elderly/adult attack ratio to 0.5-0.9 against PHIRST's 0.8. PHIRST's young/adult 1.7 (model ~1.0)
  stays unexplained. DECIDED `[owner, 2026-09-11]`: the joint model carries ONE ELDERLY
  SUSCEPTIBILITY FACTOR SHARED ACROSS COUNTRIES, with the YOUNG FIXED AT 1 -- biology is one number
  for Europe, not twelve. Prior `log2 sigma_eld ~ N(1, 0.5)` (centre 2x, 95% band 1.0-4.0x). The
  young factor is not fitted at all: the grid showed the weekly wave shapes reject it in 10 of 12
  countries, so a free young parameter would only re-absorb reporting level. Per-country fits keep
  `susc_by_age = FALSE`; `susc_fixed = c(1, 1, sigma_eld)` supplies the shared factor when the outer
  stage holds it fixed. The young/adult attack-rate gap against PHIRST is therefore a KNOWN,
  DOCUMENTED shortfall of the model, not a parameter -- report it, do not fit it (see E5).

## 4. Vaccination

- D1 `[stan][owner]` A daily vaccination fraction `delta_vax[day, a]` moves people from `S_u -> S_v` and
  `R_u -> R_v` in proportion to compartment size (infectious people are not vaccinated). Mass is
  conserved (moved amounts computed before the update; bug fixed 2026-08).
- D2 `[stan][owner]` Timing "as is": the season's entire 65+ coverage is applied as ONE pulse on
  1 October; young and medium groups receive no vaccination (no data). `[data]` Coverage: observed
  65+ history 2012/13-2021/22 (`vaccination_history_65plus`, 30 countries), post-COVID 65+ from
  `data/external/vaccination_coverage_65plus_postcovid.csv`, else the RespiCompass scenario midpoint.
- D3 `[notes]` Vaccine effects: `ve_inf = 0.25` (against infection given exposure, i.e. susceptibility),
  `ve_ili_cond_inf = 0.20` (against ILI given infection), `ve_spread = 0.20` (onward transmission).
  `[owner]` "we can do better by scraping": season-specific VE against infection from
  `data/external/vaccine_effectiveness.csv` (via `analysis_helpers.R::ve_vs_dominant()`) is a
  settings switch; default = the fixed notes values so the model matches the Stan runs.
- D4 `[stan]` Vaccinated compartments reset to zero at every season start (annual vaccination, no
  carry-over of vaccine immunity).
- D5 `[stan]` Vaccination reduces ILI given infection for the vaccinated only through
  `ve_ili_cond_inf`; it does not change the reporting proportion otherwise.

## 5. Seasons and initial conditions

- E1 `[stan]` Each season (1 Aug - 31 Jul) is an independent epidemic: at the season start the
  compartments are RESET to the initial state; natural immunity does not carry over between seasons
  (the previous season's attack rate enters only through whatever the fitted `S0` absorbs).
- E2 `[owner]` What varies where: the initial susceptibility `S0_c` is per COUNTRY and SHARED ACROSS
  SEASONS (and across the three age groups: initial R fraction = 1 - S0_c - I0 in every group; an
  age-specific modifier is a later extension). Between-season variation within a country is carried
  by the season transmissibility `R0_s` (shared across countries), the vaccination schedule, and the
  process noise -- NOT by a per-season susceptibility. (The per-season susceptibility remains the
  job of the single-population susceptibility methods; this model asks how much of the season
  variation a Europe-wide transmissibility factor explains.) The Stan model as committed had ONE
  initial simplex for all seasons, ages and -- per fit -- one country; this design keeps the shared
  initial state per country and adds the season factor on transmission instead of the disabled
  `i_season`/`r_season` initial-state deviations.
- E3 `[proposal, revised 2026-09-10 on evidence; prior recentred by the owner]` Seed:
  `I0_s = 10^-6.5 * exp(delta_s)` of each age group on day 1 of the season, with a PER-SEASON deviation
  `delta_s ~ N(0, 3)` (settings `I0_fraction`, `I0_by_season`, `prior_logI0_sd`). The centre was 1e-5
  with sd 2.5 until the five-country fits (35 seasons) put 29% of the fitted seeds below that prior's
  lower 2.5% bound (median 10^-6.4, range 10^-8.8 .. 10^-3.6): LATE waves need SMALL seeds, and a
  prior that excludes them delays every late season. The 95% band of the new prior, 10^-9.1 .. 10^-3.9,
  covers every fitted seed so far. Why not the fixed seed first agreed: with `S0_c` shared across seasons (E2) the
  fixed 1 August seed left the model no handle on WHEN a season arrives, and the first Danish fit
  showed it -- the deterministic SIR peaked around week 15-22 in every season against observed peaks
  at week 30-35, and the filter could only follow the data by abandoning the SIR (q = 0.32, phi = 1).
  The observed rise steepness matched the fitted growth rate, so the failure was timing, not
  transmission. A seed k times smaller arrives log(k)/r days later at growth rate r and leaves
  r = gamma*(R0_s*S0_c - 1), hence the identification of R0_s and S0_c, untouched -- this is the
  smooth form of the Stan model's disabled `i_season` term (prior `logit(I0) ~ N(logit 3e-6, 5)` in its
  comments). The single-population susceptibility methods never met this problem because their
  per-season S0 absorbed timing. Set `I0_by_season = FALSE` to reproduce the fixed-seed variant.
- E4 `[data]` The four COVID seasons are excluded as OUTCOMES (same rule as the panel); the season
  window and the RespiCompass/ERVISS stitch follow `stitch_iliplus.R`.
- E5 `[owner, 2026-09-11]` AGE-SPECIFIC INITIAL IMMUNITY (an `S0_a` modifier) is PARKED -- deliberately
  not fitted, to keep the model from over-complicating and to avoid a third age parameter that trades
  off against the two we keep (reporting offsets F2 and the shared elderly susceptibility C7). It is
  parked, NOT rejected: it remains the leading mechanistic candidate for the one thing the model does
  not reproduce, PHIRST's young/adult infection ratio of ~1.7 against the model's ~1.0 (C7). The
  consequence is stated as a limitation: `S0` is shared across age groups, so the model's ABSOLUTE
  attack rate in the young is probably too low, and age-group attack rates should be read as relative
  within an age group across seasons -- exactly the same caution the single-population methods carry
  for `S0` across seasons. Revisit only if a question needs the young's absolute burden.

## 6. Observation model

- F1 `[stan][data]` Observed quantity: weekly ILI+ per age group = age-specific ILI consultation rate x
  influenza positivity. `[owner]` Detection is AGE-INVARIANT: positivity is only reported for all ages
  and is applied to every band, and the reporting proportion (F2) carries no age offset (the Stan
  model's `prop_ili_age` offsets are dropped by decision). Consequently age differences in observed
  ILI+ must be explained by the DYNAMICS (contact structure, susceptibility) -- a strong, testable
  assumption.
- F2 `[stan][owner, confirmed 2026-09-11]` Expected ILI+ per age and week:
  `mu[t,a] = c_{c,s,a} * (new infections in week t, age a, vaccinated infections weighted
  (1 - ve_ili_cond_inf)) + b_source * N_a / 1e5`, with
  `c_{c,s,a} = c_c * 2^(off_a) * exp(delta_s)`: a country reporting proportion (the Stan `prop_ili`),
  AGE OFFSETS (the Stan `prop_ili_age`; medium = reference; `off ~ N(0, 1)`; settings `c_by_age`) and
  a PER-SEASON DEVIATION (the disabled Stan `prop_ili_season`; `delta ~ N(0, 0.5)`; `c_by_season`).
  Both switches default to ON on Danish evidence and await the owner's decision: (i) with age-invariant
  reporting the medium group is over-predicted ~2x in every season and the age-offset fit is 118 nats
  better (fitted: young 2.1x, elderly 2.8x the medium rate per infection); (ii) with S0 shared and R0_s
  the only season factor the 2015/16, 2017/18 and 2024/25 peaks (2-4x larger) cannot be reproduced --
  a larger R0 makes a wave sharper and earlier as well as bigger -- and the season deviation absorbs
  them (89 nats; deviations 0.5-1.7x). CONFIRMED by the owner 2026-09-11: the AGE OFFSETS are FITTED,
  one pair per COUNTRY, SHARED ACROSS that country's SEASONS -- reporting is a property of a
  surveillance system, not of a season. DECIDED `[owner, 2026-09-11]`: the SEASON DEVIATION
  `exp(delta_s)` is ONE VALUE PER SEASON SHARED ACROSS COUNTRIES in the joint model (8 parameters,
  not 96) -- the observation-side twin of the shared `R0_s`: strain symptomaticity is viral and
  Europe-wide, and the per-country deviations already moved together (FR-ES 0.97, FR-NO 0.82,
  DK-ES 0.78). Interpretation: ILI+ per infection in season s relative to the norm, a mix of strain
  symptomaticity (shared), care-seeking/testing that year (local, now absorbed elsewhere or into
  phi) and season size the mechanism cannot produce. Country-season misfit therefore lands in phi,
  S0_c or R0_s and must be watched (Estonia is the known outlier). Per-country fits keep it free
  per season as the warm start and the sharing diagnostic. The baseline `b` is PER SOURCE (`b_by_source`): RespiCompass
  ILI+ is exactly zero in weeks without detections (24-52 zeros per pre-COVID season in DK) while the
  ERVISS reconstruction has a positive floor, and one shared b forced phi towards 1.
- F3 `[stan][proposal]` Scale: ILI+ rates per 100 000 of the age group are converted to COUNTS via
  `rate * pop_a / 1e5` -- with EACH age group's own population (the legacy `make_stan_list()` indexed
  `pop_age_group[1,]`, i.e. scaled every age group by the 0-4 population; a bug, not replicated) --
  so that the negative-binomial-like variance `Var = mu + mu^2/phi` is on a count scale; the Kalman
  innovation uses this variance with a Gaussian innovation likelihood (the EKF's prediction-error
  decomposition; the exact NB2 of the Stan model is not available inside a Kalman filter). The
  conversion is a SCALE device: the reporting proportion `c` absorbs it, so countries whose ILI rate
  is per 100 consultations (CY, LU, MT, x1000 in the panel) or per 100 000 consultations (FI) are
  usable, only their count-scale variance is nominal. Missing weeks are skipped (no update), not
  zero-filled. Counts are not rounded (the innovation is Gaussian).
- F4 `[owner]` NO likelihood weights (`weight_obs_epi`) and NO separate cumulative-burden likelihood
  term (`n_season_cum_fit`, `sigma_cum_ili`): the weekly likelihood already contains the season's
  cumulative burden, and down-weighting the weekly points by 0.1 to balance a duplicated term
  tempered the likelihood (inflating uncertainty ~sqrt(10)) without a principled basis.
- F5 `[data]` RespiCompass (<= 2023/24) and ERVISS (>= 2023/24) age-specific ILI+ are on one scale via
  the per-country alignment factor of the panel stitch (estimated on the age TOTALS and applied to
  every band); per-100-consultations countries (CY, LU, MT) are scaled x1000. Country-seasons enter
  if and only if they are in the committed panel (the total-based inclusion rule); a band with
  missing weeks is simply not updated in those weeks. `[data]` Age-complete series (>= 15 positive
  weeks in all four bands) exist for all 8 panel seasons in DK, EE, ES, FR, NO and for >= 6 seasons
  in 12 countries; HU, GR and LU report age totals only in both sources and cannot be age-structured.
- F6 `[data]` Vaccination coverage for the young and medium groups is assumed ZERO. A sparse
  non-elderly coverage file exists (`data/vax_flu_history_all.csv`, loaded as
  `data$vax$data_vax_history_all` but consumed by nothing: an Italian age ladder, NO/SE 0-17 and
  18-64, 18+ for nine countries in 2018/19-2020/21); it is a later option, not a default.
- F7 `[docs]` `documentation/documentation.Rmd` diverges from the Stan code in three places and this
  model follows the CODE: (i) it names the vaccine effects differently (its `ve_inf` = infectiousness
  = code `ve_spread`; its `ve_susc` = code `ve_inf`; its `ve_severe` = code `ve_ili_cond_inf`);
  (ii) its vaccination step moves `vax * S_u/(S_u+R_u)` people, the code moves the fraction
  `delta_vax` of each of S_u and R_u; (iii) its I_v equation lacks the `(1 - ve_susc)` factor. The
  Rmd is to be brought in line when this model is documented.

## 7. Process noise (the Kalman filter)

- G1 `[owner]` The filter's purpose is to let the fitted "true" trajectory wiggle more than a
  deterministic SIR fitted to the data would -- ordinary state process noise, NOT an autocorrelated
  transmission process.
- G2 `[proposal, revised on evidence]` Process noise enters the infected compartments
  MULTIPLICATIVELY (sd = q x I per week), FIXED at q = 0.05 (settings `q_fixed`; NULL estimates it under
  `log q ~ N(log 0.1, 0.5)`). Why fixed: estimated freely on Denmark the filter escaped to q = 0.3-0.6,
  at which point it carried the wave by state corrections and the mechanistic parameters stopped
  mattering to the likelihood (the 'loose filter masks the model' failure in decisions.md). In an
  exponentially growing system even a small weekly q compounds, so q is a modelling choice ('how much
  may the trajectory wiggle'), not a well-identified parameter. Proportional noise is scale-free
  across age groups and seasons and avoids the early-rise looseness of additive noise.
- G3 `[proposal, revised on evidence]` NO initial-state uncertainty (P0 = 0). The single-population
  EKF used sd = 5% of S0 and I0; here that covariance scales with the fitted per-season seed, and the
  optimiser exploited it: the EKF objective genuinely preferred a degenerate mode (negll 6564 vs 9100
  at the deterministic optimum) in which an inflated seed inflated P0 and the filter did the fitting
  while c, b and phi collapsed. With P0 = 0 the ordering reverses (5573 vs 7366) and the EKF optimum
  agrees with the deterministic stage in every parameter (S0 0.830 vs 0.834; R0_s within 0.02).
- G4 `[proposal]` TWO-STAGE fitting (settings `two_stage`): stage 1 fits the deterministic model
  (y ~ N(mu, mu + mu^2/phi), no filter, multi-start); stage 2 starts the EKF there. Stage 1 is the
  well-posed objective that pins timing and shape; stage 2 is the noise-aware refinement. Both
  optima are kept in the fit object and drawn in the fit figure (dotted = stage 1).
- G5 `[proposal]` After each Kalman update the state is clamped to [0, 1] with the floor at EXACTLY 0,
  as in the deterministic step. (An earlier positive floor of 1e-12 re-seeded every compartment every
  week -- the vaccinated infectious before the pulse in particular -- and made the filter drift off the
  deterministic model by ~1e-12/I0 per week; noticeable once seeds of 10^-8 are admitted. The core
  test now requires the zero-noise EKF to reproduce the deterministic model to 1e-9 at a 1e-9 seed.)

## 8. Priors (as optim penalties on transformed parameters)

- H1 `[proposal]` `logit(S0_c) ~ N(logit 0.75, 1)` (from the susceptibility methods).
- H2 `[owner][proposal]` `log(R0_s) ~ N(log 1.5, 0.05)`. Sensitivity on Denmark: widening to sd 0.1-0.3
  changed nothing (the likelihood itself keeps R0_s within 1.44-1.60), so the width is not what limits
  the fit; season SIZE differences are carried by the reporting deviation (F2), not by R0.
- H3 `[proposal, revised on evidence]` `log(phi) ~ N(log 4, 0.5)`. The age-specific weekly ILI+
  series scatter 23% (median) around a 3-week moving average, i.e. phi ~ 3 is the irreducible noise
  even for a perfect mean; the earlier centre of 15 was wrong for these data. `log c ~ N(log 0.05, 2)`
  (weak; closes a 'no epidemic' trap at c -> 0); `b` unpenalised.

## 9. Inference

- I1 `[owner]` Staged: (i) a base-R REFERENCE implementation (readable, tested against the
  deterministic simulator, against a growth-rate check of the R0 calibration, and against synthetic
  parameter recovery); (ii) per-country fits: `S0_c` + regularised per-season `R0_{c,s}` + c, b, phi, q,
  multi-start optim on the EKF likelihood (the existing harness); (iii) an Rcpp port of the Kalman
  loop, and (iv) the JOINT fit of all countries and seasons with `R0_s` shared (feasibility note in
  the settings file).
- I2 `[owner]` The base-R model and the C++ implementation must remain IDENTICAL:
  `tests/testthat/test-comp-model-cpp.R` evaluates both engines on the same data and parameters
  (synthetic with missing cells, a pulse and large process noise; the real Danish fit when present)
  and requires the log-likelihoods, one-step-ahead means and filtered S and I to agree to 1e-10
  relative. `comp_model_core.cpp` carries the header note that the R version
  (`comp_model_core.R`) is the reference and that every model change must be made in both files.
  Measured speed-up ~17x; the fitter and the joint fit dispatch through
  `cm_ekf_season_engine(engine = "R" | "cpp")` (`comp_model_cpp.R`).
- I3 `[proposal]` Point estimates are MAP (penalised EKF likelihood); uncertainty from the Hessian
  (Laplace) at the optimum; no MCMC.

## 10. Deliberately NOT modelled

Waning, ageing, births/deaths, importation after the seed, latent period, multiple co-circulating
strains (documented single-strain caveat), age-specific positivity, vaccination of < 65 year-olds,
NPIs, COVID seasons as outcomes, care-seeking changes between eras (absorbed by `c` per country and
the source alignment).

## 11. Open decisions (`[open]`)

1. ~~The AGE MECHANISM~~ -- SETTLED `[owner, 2026-09-11]`: reporting offsets per country (F2) plus
   ONE elderly susceptibility factor shared across countries, young fixed at 1 (C7). Age-specific
   initial immunity is PARKED, not rejected (E5).
2. ~~Per-season drift of the reporting proportion `c`~~ -- SETTLED `[owner, 2026-09-11]`: one season
   deviation per season, SHARED ACROSS COUNTRIES in the joint model (F2).
3. Prior widths in H1-H3 and the process-noise prior (G2).
4. Season-specific VE from the CSV vs the fixed notes values (D3).

## 12. The parameters in words (a glossary for reading the fit figures)

- `c` (F2) -- the REPORTING proportion: expected ILI+ per infection, i.e. P(ILI consultation | infection)
  x P(sampled and tested) x positivity, on the ILI+ rate scale of the panel. Per country; multiplied by
  the age offsets `2^off_a` (young, elderly vs the medium reference) and the season deviation
  `exp(delta_s)`. It is what turns modelled infections into the observed series, so a country with
  low `c` may have more infections than a country with a higher curve.
- `b` (F2) -- the off-season BASELINE of the observed series: ILI+ per 100 000 that is NOT influenza
  infections in the model's sense (residual positivity in ILI from other causes, misclassification,
  the reconstruction floor of the stitched series). There is one `b` per DATA SOURCE (`b_RespiCompass`,
  `b_ERVISS`) because the two sources have different floors: RespiCompass ILI+ is exactly zero in
  weeks without detections, the ERVISS reconstruction sits on a small positive floor -- one shared
  `b` cannot be right for both and the mismatch was absorbed by phi.
- `phi` (F3) -- OBSERVATION (measurement) noise: the overdispersion of the weekly ILI+ around the
  model's expected value, `Var = mu + mu^2/phi`, so 1/sqrt(phi) is the coefficient of variation at
  large counts (phi = 4 means +/- 50%). Mechanistically it is everything that scrambles the
  MEASUREMENT of a given week's infections without changing the epidemic: sentinel sampling, the
  number of swabs tested, holiday reporting, weekly positivity sampling error in a small
  denominator, age-band reallocation between sources. It does NOT propagate: a noisy week leaves no
  trace on the next week's state.
- `q` (G2) -- PROCESS noise: the weekly wiggle of the latent epidemic itself, a multiplicative sd of
  `q x I` on the infected fractions (fixed at 5%). Mechanistically it stands for what the SIR leaves
  out of the TRANSMISSION process: school holidays and weather changing contacts for a week,
  importations and spatial heterogeneity (a second city taking off), a co-circulating strain, the
  stochasticity of a small number of chains early in the season. Unlike phi it PROPAGATES: a
  perturbation of I this week is carried forward by the dynamics and changes every later week -- which
  is why the filter must be kept tight (G2-G3): a loose q lets the filter carry the wave by state
  corrections and the mechanistic parameters stop mattering.
- `S0` (E2) -- the fraction of each age group susceptible at the season start (per country, shared
  across seasons and ages); `R0_s` (C4) -- the season's transmissibility at full susceptibility given
  the age mixing; `I0_s` (E3) -- the season's seed, i.e. its ARRIVAL time.
