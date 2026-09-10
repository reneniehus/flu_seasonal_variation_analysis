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
- E3 `[proposal]` Seed: `I0 = 1e-5` of each age group, fixed, planted on day 1 of the season (as the
  susceptibility methods do; `decisions.md`). With I0 and the seed day fixed, onset timing is carried
  entirely by the growth rate, which is what identifies `S0`. (The Stan model fitted a shared initial
  I fraction instead; prior `logit(I0) ~ N(logit 3e-6, 5)` in its comments.)
- E4 `[data]` The four COVID seasons are excluded as OUTCOMES (same rule as the panel); the season
  window and the RespiCompass/ERVISS stitch follow `stitch_iliplus.R`.

## 6. Observation model

- F1 `[stan][data]` Observed quantity: weekly ILI+ per age group = age-specific ILI consultation rate x
  influenza positivity. `[owner]` Detection is AGE-INVARIANT: positivity is only reported for all ages
  and is applied to every band, and the reporting proportion (F2) carries no age offset (the Stan
  model's `prop_ili_age` offsets are dropped by decision). Consequently age differences in observed
  ILI+ must be explained by the DYNAMICS (contact structure, susceptibility) -- a strong, testable
  assumption.
- F2 `[stan][proposal]` Expected ILI+ per age and week: `mu[t,a] = c_country * (new infections in
  week t, age a, with vaccinated infections weighted (1 - ve_ili_cond_inf)) + b_country`, with `c`
  the reporting proportion (per country, shared across seasons and ages; the Stan `prop_ili`) and `b`
  an off-season baseline (from the susceptibility methods; the Stan model had none).
  `[open]` whether `c` may drift by season (a weak per-season deviation, as the disabled Stan
  `prop_ili_season` intended).
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
- G2 `[proposal]` Process noise enters the infected compartments MULTIPLICATIVELY (on log I, i.e. an sd
  proportional to the current I) with one shared variance `q` per fit, regularised small
  (`log q ~ N(log 0.05, 1)` `[open]`). Rationale: additive noise on I is loosest, in relative terms,
  exactly during the seeding/rise window that identifies `S0` (2026-08 review); proportional noise is
  scale-free across age groups and seasons.
- G3 `[stan-r]` Initial-state covariance is tight (as in the existing EKF: sd = 5% of S0 and I0), so
  the fitted initial condition is trusted and the noise only makes modest weekly corrections.

## 8. Priors (as optim penalties on transformed parameters)

- H1 `[proposal]` `logit(S0[c,s]) ~ N(logit 0.75, 1)` (from the susceptibility methods).
- H2 `[owner][proposal]` `log(R0_s) ~ N(log 1.5, 0.05)`.
- H3 `[proposal]` `log(phi) ~ N(log 15, 0.8)`; `log(q)` as G2; `c`, `b` unpenalised (data-scale).

## 9. Inference

- I1 `[owner]` Staged: (i) a base-R REFERENCE implementation (readable, tested against the
  deterministic simulator, against a growth-rate check of the R0 calibration, and against synthetic
  parameter recovery); (ii) per-country fits: `S0_c` + regularised per-season `R0_{c,s}` + c, b, phi, q,
  multi-start optim on the EKF likelihood (the existing harness); (iii) an Rcpp port of the Kalman
  loop, and (iv) the JOINT fit of all countries and seasons with `R0_s` shared (feasibility note in
  the settings file).
- I2 `[owner]` The base-R model and the C++ implementation must remain IDENTICAL: a test evaluates
  both likelihoods and both filtered trajectories on the same data and parameters and requires
  agreement to numerical precision (1e-10 relative); the C++ code carries a header note that the
  R version is the reference and every change must be made in both. (Test and note are added with
  the port, stage iii.)
- I3 `[proposal]` Point estimates are MAP (penalised EKF likelihood); uncertainty from the Hessian
  (Laplace) at the optimum; no MCMC.

## 10. Deliberately NOT modelled

Waning, ageing, births/deaths, importation after the seed, latent period, multiple co-circulating
strains (documented single-strain caveat), age-specific positivity, vaccination of < 65 year-olds,
NPIs, COVID seasons as outcomes, care-seeking changes between eras (absorbed by `c` per country and
the source alignment).

## 11. Open decisions (`[open]`)

1. Age-specific `S0` modifier (E2) -- start shared, revisit after the first fits.
2. Per-season drift of the reporting proportion `c` (F2).
3. Prior widths in H1-H3 and the process-noise prior (G2).
4. Season-specific VE from the CSV vs the fixed notes values (D3).
