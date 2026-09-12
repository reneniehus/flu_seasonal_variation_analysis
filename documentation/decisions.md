# Design decisions & rationale

The README records *what* the repo does and how to run it; this file records *why* the key
modelling, method and data choices were made, so the reasoning survives beyond commit messages and
chat. Append new decisions as they are made. Keep each entry short: the **decision**, the reason,
and the main alternative considered.

## Inference & modelling

- **Fix R0 from the literature (`susc_R0 = 1.5`); fit a single susceptibility `S0` per season.** How
  fast a season's virus spreads is set by the effective reproduction number `R_eff = R0 * S0`
  (equivalently the early growth rate `r = gamma*(R0*S0 - 1)`): intrinsic transmissibility `R0` and
  the susceptible fraction `S0` are multiplicatively linked, so a single curve identifies their
  *product*, not each separately — fitting both jointly is degenerate (R0 drifts to its prior and S0
  follows). We therefore fix R0 and let one per-season `S0` act as the "sensor" of how easily that
  season's virus spread. This is appropriate because the factors of interest here act mainly on the
  *susceptible fraction* — vaccine coverage, antigenic match/mismatch and prior-season population
  immunity (see `PROJECT_SCOPE.md`, refs 4, 6–8) — whereas the determinants of intrinsic
  transmissibility (contact patterns, behavioural restrictions, viral subtype infectivity) are either
  out of scope or reasonably treated as season-invariant. Inferring a susceptibility state from
  incidence under a fixed transmission model is the established *susceptible-reconstruction* approach
  (Bjørnstad et al. 2002; Finkenstädt & Grenfell 2000). **Caveat:** holding R0 fixed means any genuine
  season-to-season variation in transmissibility — notably subtype infectivity, with A(H3N2) seasons
  more transmissible than A(H1N1) (Infection, Genetics and Evolution 2018) — is absorbed into the
  fitted S0. The fitted S0 is therefore a *composite* index of "how easily the season spread", not a
  pure immunity measure; in particular, an apparent association between S0 and dominant subtype may
  partly reflect this absorption rather than susceptibility alone.

- **Fix the seed `I0` (`susc_seed_i0 = 1e-5`), plant it at the season start (~Aug), and read the
  reporting fraction `c` only in relative terms.** The seed size `I0` and the reporting fraction `c`
  both act on the *scale* of the modelled curve — a larger seed and a larger reporting fraction are
  largely interchangeable in their effect on the observed level — so they are not jointly identifiable
  from one curve. We resolve this by fixing `I0` (a small constant import, encoding the assumption
  that flu is reliably re-seeded each year, e.g. from the southern hemisphere) and fitting `c`.
  Planting the seed at a *fixed time* (the season start) as well as a fixed size means a wave's timing
  is then explained by `S0` alone, so "an earlier / steeper season = a more susceptible population"
  holds. Consequence: the *absolute* value of `c` is not interpretable (it is anchored to the
  arbitrary `I0`), whereas *relative* variation in `c` across seasons may still carry meaning
  (care-seeking, testing intensity, severity-driven detection). *Alternative:* fit `I0` per season —
  rejected (confounds the curve's scale and timing with `c` and `S0`).

- **Interpret S0 as RELATIVE across seasons, not absolute.** The absolute level is conditional on the
  fixed R0 and seed (e.g. a lower R0 forces a higher S0); only the ranking / spacing across seasons
  is robust — and that is the quantity of interest.

- **Per season fit S0 and the reporting fraction `c`; share the baseline `b`, overdispersion `phi`
  and (EKF) process noise `q_I` across a country's seasons.** S0 (a rate) and c (a level) control
  different features of the curve, so they separate cleanly, and both are allowed to vary season to
  season (susceptibility; reporting via severity / testing / behaviour). `b`, `phi` and `q_I` are
  instead **shared across the country's seasons**: they are treated as stable country / reporting /
  process properties, and sharing them keeps the parameter count down and **stabilises the fit**.
  *Alternative:* a single shared c (over-constrains size vs timing), or per-season `b` / `phi` / `q_I`
  (more parameters competing with S0 and c).

- **Off-season activity is an additive observation baseline `b`, not part of the dynamics.** The
  off-season ILI+ floor is non-SIR (sporadic / cross-reacting detections); absorbing it in `b` frees
  the SIR to fit the wave. *Alternative:* restrict the likelihood to the epidemic window — rejected
  as circular (the window depends on the fit we are trying to make).

- **`onset_week` is read off the fitted curve with a simple threshold, and kept as a model feature.**
  For the EKF the curve is the filtered mean, which tracks data and can give early onsets; this is
  accepted as a property of that method rather than something to "fix".

- **Neg-binomial-like observation noise** (`Var = mu + mu^2/phi`). Inherited from the Stan generative
  model — but note that Stan fits absolute *counts* whereas we fit a *rate* (ILI+), so a count-style
  variance may not be ideal here. Alternative noise models (e.g. multiplicative / log-scale) should be
  tested — see the planned sensitivity analyses below.

- **Fix the recovery rate `gamma` from a 3-day mean infectious period** (`susc_infectious_period_days`).
  A defensible literature value, fixed for simplicity; estimates span ~2–4 days, and `gamma`
  propagates into `beta = R0*gamma` and the rise-rate scaling, so the value moves the fitted numbers.
  Kept fixed, but flagged for the planned sensitivity analyses below.

## Method framework

- **A swappable multi-method framework** (`sir_core.R` + `methods/` + `methods_registry.R`): every
  method fits the same per-season panel and returns one common summary table — the input the
  downstream correlation analysis consumes. Adding a method = one new file + one registry line.

- **Primary-analysis method choice is deliberately deferred.** The project is still maturing, so the
  entries below record the *rationale and observed properties* of each method, not a final selection
  of which method(s) lead the primary analysis. All methods are kept and compared; the
  phenomenological vs mechanistic comparison is itself part of the inquiry.

- **Properties of the deterministic SIR (the rigid end of the spectrum).** The deterministic method
  treats the SIR as the truth with residuals as pure noise, which makes it the most rigid: it cannot
  bend to seasons whose shape departs from a single SIR (e.g. the early, sharp FR 2022/23 wave fails
  outright). The **EKF** sits at the other end, admitting genuine process noise (a fitted `q_I` that
  lets the latent state wander off the SIR). Which of these — if either — is carried into the primary
  analysis is left open (see the deferral note); at minimum the deterministic fit is a transparent
  reference / sanity check.

- **EKF process noise is fitted but regularised small, with a tight initial covariance.** The aim is
  to keep S0 identified by the rise rate rather than absorbed by filter "tracking": a small initial
  covariance `p0` and small process noise `q_I` keep the latent-state covariance small relative to the
  observation noise, so the filter trusts the dynamics over the data. **Caveat (finite sample):** with
  only ~50 weekly observations per season the state covariance still accumulates enough over the
  season that S0 is not fully pinned — seen as the EKF drawing more between-season S0 contrast
  (0.70–0.90) than the deterministic fit (0.73–0.76), while the descriptive rise-rate agrees with the
  deterministic. So read the EKF *ranking*, not the absolute gaps. **This limited-data caveat must be
  stated clearly in the final reporting**, and — if the EKF stays in focus — probed by a sensitivity
  analysis on `p0` / `q_I` (see the planned sensitivity analyses).

- **Take the descriptive (phenomenological) characterisation as a first-class, low-assumption lens —
  and do not over-trust the mechanistic SIR.** A single-population SIR (deterministic or filtered)
  assumes the national ILI+ curve is generated by *one well-mixed epidemic*. National surveillance
  curves are not: they are the **overlay of many local epidemics** that ignite at different times and
  travel across a country, and the aggregate signal is further blurred by reporting delays that
  themselves vary in space and time. Influenza spreads hierarchically, as travelling waves structured
  by human mobility, taking weeks to months to diffuse nationally (Viboud et al. 2006; Gog et al.
  2014). A mechanistic SIR fitted to the aggregate therefore **mis-specifies the data-generating
  process**, and its fitted parameters (e.g. `S0`) can be distorted by aggregation — inviting spurious
  mechanistic interpretation. Following the precautionary, parsimonious principle, we judge it more
  honest to **compare the observable shapes** of the epidemics directly (onset, peak, intensity,
  steepness, burden) and look for predictable patterns than to assume a transmission mechanism we know
  holds only approximately. For an EU agency working closely with the member-state experts who produce
  these data, this is also the more **responsible** stance: phenomenological summaries make minimal,
  transparent assumptions about each country's data instead of imposing a model that may not hold
  uniformly across settings. We therefore treat the phenomenological characterisation as a first-class
  approach alongside the mechanistic fits, and keep all methods; which method(s) lead the primary
  analysis is left open at this stage (see the deferral note), and comparing the phenomenological and
  mechanistic readings is itself part of the inquiry.
  *Method specifics:* smooth each curve and extract AUC, peak height, onset and steepness — all
  observed-shape descriptors. None is a clean susceptibility: the observed national rise reflects not
  only transmissibility but also how a country's local epidemics overlay (staggered starts, travel)
  and how reporting is delayed in space and time. AUC, peak height and steepness are therefore
  compared across seasons *within* a country, not across countries (onset timing likewise needs care
  across countries). The method deliberately does **not** map steepness onto an SIR susceptibility —
  that mechanistic-vs-phenomenological mapping is a separate question, out of scope here.

- **Smoothing = centered moving average, window 4 (not loess).** On regular, single-wave weekly data
  a moving average gives essentially the same features as loess (steepness ~0.31 vs 0.29), so the
  simpler, more transparent option wins. The even "2×4" window (half-weights at the ends → no timing
  shift) preserves the peak better than window 5 while staying smooth. *Alternative:* loess (fine but
  less parsimonious) or a GCV spline (numerically fragile on these short series).

## Data & infrastructure

- **Four seasons impacted by the acute COVID-19 pandemic phase are excluded:** 2019/2020, 2020/2021,
  2021/2022 and 2022/2023. They are treated as disrupted by the pandemic and the response to it (NPIs,
  collapsed/atypical influenza circulation, and changes in care-seeking and testing). The analysis
  therefore spans 2014/2015–2018/2019 and 2023/2024–2025/2026 (8 seasons).

- **The two ILI+ sources span different eras and must be combined for a long record.** ERVISS sentinel
  ILI+ (reconstructed as ILI consultation rate × influenza positivity) only reaches back to 2020/21,
  whereas RespiCompass ILI+ runs 2014/15–2023/24. Covering the pre-pandemic seasons (2017/18, 2018/19)
  *and* the most recent ones (2024/25, 2025/26) therefore requires both sources, with **2023/24 as the
  overlap** for cross-checking. Because the two construct ILI+ differently, mixing them along the time
  axis risks a spurious "source/era effect" in the very inter-seasonal comparison of interest; the
  overlap season is used to quantify and adjust for that difference. (See `code/03_report/
  data_availability.R` for the coverage picture.)

- **ERVISS ILI+ is reconstructed RespiCompass-style.** ILI+ = ILI consultation rate × influenza test
  positivity, using **sentinel** positivity (re-derived as detections/tests), **except non-sentinel
  positivity for Malta, Iceland, Croatia, Romania, Latvia and Finland**, which lacked an adequate
  number of weeks of sentinel test positivity — following RespiCompass. This assumes non-sentinel test
  positivity reflects influenza positivity in a similar way as sentinel data would. **Consultation-rate
  units also differ by country:** ERVISS reports per 100 000 population, except Cyprus, Luxembourg and
  Malta (per 100 consultations) and Finland (per 100 000 consultations); to match RespiCompass the
  per-100-consultations countries (CY, LU, MT) are scaled ×1000 onto the per-100 000 basis.

- **The 2023/24 overlap validates the reconstruction, explains the offsets, and anchors the stitch.**
  On 2023/24, **15 of 24 countries match RespiCompass exactly** (cor ≈ 1, ratio ≈ 1). The deviations
  are understood, not method errors:
  - **Consultation-rate units (CY, LU, MT):** off by ×1000 (per-100-consultations vs per-100 000) — a
    deterministic, principled correction; after it LU matches exactly.
  - **Positivity-construction differences (IE ≈ ×0.59, BE ≈ ×1.18, LT):** RespiCompass's influenza
    positivity differs from the re-derived detections/tests by a roughly constant per-country factor;
    once applied the dynamics are identical (IE cor 1.00, BE 0.996; LT 0.91, a few noisy weeks).
  - **LV → ERVISS only (decided).** The RespiCompass LV series is anomalous (flat all winter, a single
    May 2024 spike) while the ERVISS reconstruction is sensible (a January peak) — the source data, not
    the method, is wrong. LV's RespiCompass seasons are therefore dropped and LV is kept on the ERVISS
    reconstruction only (2023/24–2025/26), treated like the single-source countries below.
  - **No clean overlap (NO, ES, SK):** kept on a **single source — whichever has the most non-COVID
    seasons — dismissing the other**: NO → RespiCompass (6 seasons), ES → RespiCompass (5), SK → ERVISS
    (3). No alignment factor is applied to single-source countries (LV likewise: ERVISS only).
  **Stitch:** RespiCompass ≤ 2023/24 + ERVISS reconstruction 2024/25+, with a per-country alignment
  factor from 2023/24 applied to the ERVISS era (= 1 for the 15 exact-match countries; the unit /
  positivity factors above otherwise).

- **A committed slim panel (`data/slim_flu_iliplus.csv`), loadable in base R.** The susceptibility
  fits run offline from this file (a contiguous weekly grid, seeded from the season start), so they
  need no tidyverse pipeline — fast, dependency-light, reproducible. Countries: DK, FR, IE, HU (each
  complete across all 8 non-COVID seasons; chosen for completeness), set via `params$susc_countries`.

- **Retired the single-season, free-R0 EKF tracking fit.** Superseded by the EKF susceptibility
  method; removing it keeps the repo focused. The shared `.ekf_filter` engine is kept in `sir_core.R`.

- **Retired the legacy SIR run-path scaffolding.** The original-pipeline orchestration was dead and
  fully superseded by the method framework (`sir_core.R` + `methods/` + `methods_registry.R`) for the
  driver analysis. Removed: the no-op `run_model.R` stub, the unreachable `process_and_save.R`
  (RespiCompass hub-submission saver), the five orphaned single-model scripts (`model_SIR_multiseason`,
  `model_SIR_simple`, `model_SIR_simple_r0`, `model_arima_simple`, `model_last_year_burden` — several
  also broken: missing Stan files, a stray `browser()`, a `load()`-for-`log()` typo), and the dead
  helpers they alone used in `flu_functions.R` (`data_into_all_season`, `get_contact_matrix`, the
  severity helpers, `rep_warning_wed`, `mnaming`, `squash_axis`). None was sourced or called by the
  active path (tests, `build_slim_panel.R`, `code/05_analysis/`); the reference Stan model
  `stan/SIR_multiseason_age_vax_2.stan` is kept. Git history preserves all removed code, and `00_main.R`
  is now data-build + eyeballing only.

- **NA handling.** The moving average is NA-aware and bridges short gaps; a week whose whole window
  is missing stays NA and is skipped by the feature helpers. As more countries bring longer internal
  gaps, the planned step is linear interpolation of internal gaps before smoothing (hook noted in
  `.smooth_curve`); leading / trailing off-season NA can stay (≈ 0). Not built yet — the current
  panel has no problematic gaps.

- **Environment.** CRAN is blocked in the managed environment; dependencies install from Posit
  Package Manager binaries (the setup script sets the repo URL + HTTP user agent) and are restored
  via `renv`.

- **External driver data pulled into `data/external/`.** To fill the seasons/variables missing in-repo,
  four candidate-driver streams were web-sourced (full provenance + caveats in
  `documentation/external_drivers.md`): (1) **dominant subtype** for the 5 pre-COVID seasons (ECDC /
  Eurosurveillance), completing all 8 seasons and enabling `code/05_analysis/subtype_8season.R` (subtype
  now recurs across both eras → far less season-confounded); (2) **65+ vaccination coverage** for
  2023/24–24/25 (~10 panel countries; Eurostat/ECDC/national; comparability caveats); (3) **winter
  climate** — only *continental* C3S anomalies were obtainable here (gridded reanalysis APIs are egress-
  blocked), so per-country ERA5 is a documented next step; (4) **vaccine effectiveness** (I-MOVE/VEBIS,
  point estimates + 95% CIs). Constraint of record: in this environment `WebFetch` and data APIs return
  403 (egress policy); only `WebSearch` was available, so literature data was pulled but not raw APIs.
  Decision: the dominant subtype is assumed ~identical across EU/EEA countries within a season (continental).

## Open questions & planned sensitivity analyses

Choices made for simplicity or by assumption that should be probed downstream, once the analysis
pipeline is in place:

- **Observation-noise model** — test alternatives to the neg-binomial-like variance (the data are
  rates, not counts).
- **Infectious period / `gamma`** — vary around the 3-day default.
- **Fixed `R0` value** — vary around 1.5.
- **Seed `I0` and seed week** — vary the seed size and the season-start week (the latter also tests
  the "early season = more susceptible" reading).
- **EKF `p0` / `q_I` and finite-sample identifiability** — with only ~50 weekly points per season the
  EKF S0 is not fully pinned (see the EKF entry); vary the initial covariance and process-noise
  regularisation to quantify how much of the between-season S0 spread is filter freedom rather than
  signal. To be added at the very end, if the EKF remains in focus.
- **FOLLOW-UP (surveillance/ERVISS experts): why ERVISS and RespiCompass positivity differ for IE,
  BE, LT.** Their ILI+ offset is a *positivity-construction* difference, not a units one — RespiCompass's
  influenza positivity differs from the re-derived ERVISS detections/tests by a roughly constant
  per-country factor (IE ≈ ×0.59, BE ≈ ×1.18; LT matches but is noisy). The exact reason (different
  test denominator, a different/blended positivity source, published vs re-derived) is unknown and
  should be clarified with ERVISS / surveillance colleagues. For now these countries are aligned
  empirically via the 2023/24 overlap factor.

## References

Methodological and supporting literature for the decisions above. Driver references (dominant
subtype, vaccine coverage, antigenic match / effectiveness, prior immunity) are in `PROJECT_SCOPE.md`
(refs 4–8).

- Viboud C, Bjørnstad ON, Smith DL, Simonsen L, Miller MA, Grenfell BT. Synchrony, waves, and spatial
  hierarchies in the spread of influenza. *Science.* 2006;312(5772):447–451.
- Gog JR, Ballesteros S, Viboud C, et al. Spatial transmission of 2009 pandemic influenza in the US.
  *PLoS Computational Biology.* 2014;10(6):e1003635.
- Bjørnstad ON, Finkenstädt BF, Grenfell BT. Dynamics of measles epidemics: estimating scaling of
  transmission rates using a time series SIR model. *Ecological Monographs.* 2002;72(2):169–184.
- Finkenstädt BF, Grenfell BT. Time series modelling of childhood diseases: a dynamical systems
  approach. *Journal of the Royal Statistical Society: Series C (Applied Statistics).*
  2000;49(2):187–205.
- Transmissibility and severity of influenza virus by subtype. *Infection, Genetics and Evolution.*
  2018. https://www.sciencedirect.com/science/article/abs/pii/S1567134818306051

## 2026-08 full-code-review decisions (behaviour-affecting corrections + their rationale)

A structured multi-agent review (module reviewers + adversarial verification of every substantive
finding) drove the following decisions. Behaviour-preserving refactors (shared stitch helper, shared
Gibbs helpers, setup.R pruning) are not decisions and are only noted in commit messages.

- **Dominant subtype: hierarchical, type-first rule.** The old rule counted 'B (unknown)' toward B but
  discarded 'A (unknown)' — asymmetric, and consequential: unsubtyped A dominates reported A (2024/25
  EU/EEA: ~206k of ~260k A detections), so 2024/25 was called "B" in 22/30 countries while type-level A
  exceeded B in 28/30. New rule: (1) type by plurality of ALL characterised A vs ALL characterised B;
  (2) subtype by plurality among subtyped A (dominant = NA when type A wins with zero subtyped A).
  2024/25 flips to A(H1N1) (type share 0.74, H1 share 0.60 — genuine co-circulation, confidence
  'medium'); 2023/24 and 2025/26 unchanged. The ERVISS 'EU/EEA' aggregate is no longer treated as a
  country; it feeds a separate continental cross-check table. Downstream: the "B seasons carry the
  largest burden" post-COVID claim dissolved (those were H1N1 seasons), and the 8-season "every subtype
  recurs in both eras" motivation weakened to H1N1+H3N2 only (B stays essentially 2017/18).

- **covid_era is defined by SEASON, not by data source.** post = 2023/24 onward. The old
  era-from-source coding filed 20/22 countries' 2023/24 (a post-COVID season, RespiCompass-sourced)
  under "pre", mislabelling figures and diluting the post-COVID earlier-onset contrast (3.3 vs the true
  3.9 weeks). `covid_era` and `source` now travel as separate columns in descriptors.csv; models that
  need the measurement shift absorbed use season intercepts (subtype_8season) or print both groupings
  (analyse_patterns).

- **Season-level predictors get a SEASON random intercept** (subtype_8season.R, gibbs_ri2 + lme4
  crossed REs). With country intercepts only, ~20+ countries sharing one season-level value are treated
  as independent replicates and CrIs are anti-conservative. Under the honest specification most
  subtype burden contrasts no longer exclude 0; H3N2-peaks-earlier survives. The pre-COVID whisker
  model (4 seasons) keeps country-only REs but now carries the caveat in the docs.

- **VE values are derived from the provenance CSV, not hard-coded.** analysis_helpers.R::
  ve_vs_dominant() implements the explicit preference rule (primary-care all/target-group; eos point >
  interim point > range midpoint; influenza_A/any fallback). The missing 2016/17 all-target-group row
  (25.7, I-MOVE interim, ES.2017.22.7.30464) was added to the CSV; 2018/19 becomes 37.5 (the honest
  midpoint of 32-43, replacing a hand-rounded 38); 2024/25 becomes 30 (A(H1N1) end-of-season,
  replacing an untraceable "B interim" 58 that also rested on the miscounted dominance).

- **Out-of-sample test set is no longer conditioned on coverage availability.** Protection was never an
  out-of-sample predictor, yet the test filter required post-COVID 65+ coverage — shrinking the test to
  a 10-row big-reporter subset. Unconditioned test = 35 country-seasons, and with a PERSISTENCE
  comparator (predict = last season's log AUC) added, the headline reverses: persistence RMSE 0.62
  beats the full model (0.75) and the country-only baseline (0.90); the full model's within-country
  moves anti-correlate with observation (r = -0.76, 17 paired countries) and its 95% PI covers only
  74%. Recorded as the honest finding; the earlier "RMSE 0.71 -> 0.46" headline was an artifact of the
  conditioned test set plus the 2024/25 "B" label.

- **2022/23 as predictor, never as outcome.** The COVID-season exclusion is about disrupted WAVE SHAPE
  (invalid as an outcome/fit target); the autoregressive prior-burden predictor needs the realised
  prior burden regardless. The full-panel rebuild for that purpose runs through the SAME shared stitch
  (stitch_iliplus.R) with the COVID filter off, and an in-script assertion verifies it reproduces the
  committed panel's AUCs exactly.

- **Panel stitch: the per-week gap-fill within 2023/24 is now the documented rule.** The stitch was
  documented as per-season but implemented per-week; 14 country-seasons in the 2023/24 overlap mix
  RespiCompass with aligned-ERVISS gap-fill weeks (values scale-consistent via the alignment factor).
  Decision: keep the per-week behaviour (more data, consistent scale), document it (stitch_iliplus.R
  header), and label the season by its dominant early-season source. Changing to a strict per-season
  rule would alter the committed panel for those 14 seasons for no analytical gain.

- **Contact matrices: full-block sums + Czechia recovered.** The 17->4 age-band aggregation summed only
  the lower triangle of diagonal blocks, dropping half the cross-band contact-ends (15-64 within-group
  rate understated ~1.5x, e.g. BE 7.59 -> 11.51 contacts/person); merged matrices now satisfy directed-
  ends reciprocity to machine precision. 'Czechia' is mapped to the Prem workbooks' 'Czech Republic'
  sheet (it was silently dropped), and countries genuinely absent (LI, NO) are now declared at load.
  No published result consumed models_in$contacts (the Stan model is parked), so nothing shifts.

- **EKF process noise: additive-on-I is a documented axis of the planned sensitivity analysis.** The
  fixed absolute q_I (~1e-4) is ~10x the seed I0, so in relative terms the filter is loosest exactly in
  the rise window that identifies S0 (profile-likelihood check: at the regularised q_I the data barely
  constrain S0 beyond its prior). Not hot-swapped (the priors/P0 are tuned to the additive scale);
  'additive vs proportional (log-I) q_I' joins the p0/q_I sensitivity plan above.

## 2026-09 compartmental model (`flu_comp_model` branch): assumptions fixed before the rewrite

The age- and vaccination-structured multi-season model is being rewritten in base R with an EKF.
Every assumption, with provenance (owner decision / Stan model / old notes / data / proposal), is
spelled out in `code/06_comp_model/ASSUMPTIONS.md` and encoded in
`code/06_comp_model/comp_model_settings.R`. Decisions of record made here:

- **Contact-matrix scaling.** The Stan model's `beta * a_factor * rowNormalised(C)` equals
  `(beta / cbar) * C`, whose realised R0 at full susceptibility is `R0 * rho(C)/cbar` = 1.58-1.73 across
  the 28 country matrices (median 1.65) instead of the intended 1.5, varying by country. The model now
  divides the contact matrix by its dominant eigenvalue so the next-generation matrix has spectral
  radius exactly R0_s (the normalisation this log recommended in 2026-06). `a_factor` is dropped.
- **Season-level R0_s shared across countries**, strong prior around 1.5; identified only by pooling
  across countries (within one country-season it trades off with S0), so per-country stage fits keep
  R0 fixed.
- **No likelihood tempering.** The Stan runner's `weight_obs_epi = 0.1` and the separate
  cumulative-burden term are dropped: the weekly likelihood already contains the season burden.
- **Detection age-invariant**; one reporting proportion per country; three age groups 0-14/15-64/65+.
- **Vaccination timing as before**: one 65+ pulse on 1 October, none for younger groups (no data).
- Fixed values from the owner's notes: gamma = 0.2777778/day (3.6 d), ve_inf = 0.25,
  ve_ili_cond_inf = 0.20, ve_spread = 0.20; season-specific VE from the provenance CSV is a switch.

## 2026-09 compartmental model: what the Danish pilot changed (evidence in `code/06_comp_model/ASSUMPTIONS.md`)

The first fits of the age- and vaccination-structured model on Denmark, read through the eyeballing
figures, overturned four of the assumptions first agreed. Each is a settings switch so the earlier
variant stays reproducible; the defaults now follow the evidence and await the owner's confirmation.

- **Per-season seed size (arrival time).** A fixed 1 August seed with S0 shared across seasons left no
  handle on when a season arrives: the deterministic SIR peaked at week 15-22 against observed peaks
  at week 30-35, and the filter could only follow the data by abandoning the model. A per-season seed
  (log-normal deviation) shifts arrival without touching the growth rate -- the smooth form of the
  Stan model's disabled `i_season` term.
- **Age-specific reporting (Stan `prop_ili_age`).** Age-invariant reporting over-predicted the medium
  group ~2x in every season; the age-offset fit is 118 nats better (medium adults report ~0.5x the
  young and ~0.35x the elderly per infection). This reverses the "detection age-invariant" decision,
  pending the owner: it is an age effect on ILI-per-infection or care-seeking, indistinguishable here
  from age-specific susceptibility.
- **Per-season reporting deviation (Stan `prop_ili_season`).** With S0 shared and R0_s the only season
  factor, the 2-4x larger 2015/16, 2017/18 and 2024/25 peaks cannot be reproduced -- a larger R0 also
  makes a wave sharper and earlier -- and widening the R0 prior (sd 0.05 -> 0.3) changed nothing; the
  likelihood keeps R0_s within 1.44-1.60. A season deviation on reporting absorbs the size differences
  (89 nats; 0.5-1.7x). Consequence for the joint stage: season size may need a Europe-wide severity
  factor alongside R0_s.
- **Baseline per data source.** RespiCompass ILI+ is exactly zero before the wave (24-52 zeros per
  pre-COVID season in DK); the ERVISS reconstruction has a positive floor; one shared baseline forced
  phi towards 1.
- **EKF design.** Estimating q let the filter carry the wave (q -> 0.3-0.6) and the mechanistic
  parameters wandered; a seed-scaled initial covariance was exploited the same way (the EKF objective
  preferred a degenerate baseline-plus-noise mode, 6564 vs 9100). Decisions: two-stage fitting
  (deterministic first, EKF from that optimum), P0 = 0, q fixed at 5% weekly -- the filter is a wiggle
  layer, not a replacement for the model. With these, the EKF optimum agrees with the deterministic
  one in every parameter (S0 0.830 vs 0.834; R0_s within 0.02).
- **Observation noise.** The age-specific weekly series scatter 23% around a 3-week moving average
  (phi ~ 3 irreducible); the phi prior was recentred from 15 to 4.
- **C++ engine.** Verified identical to the R reference to 1e-10 on synthetic and real data
  (`tests/testthat/test-comp-model-cpp.R`, 39 assertions); the default engine for fitting.

## 2026-09 compartmental model: which age mechanism? (12 countries, PHIRST)

**Question.** Age-invariant reporting over-predicted the medium group in Denmark. Two mechanisms can
absorb that and are not separable on one country: age-specific REPORTING (ILI+ per infection differs
by age; dynamics untouched) or age-specific SUSCEPTIBILITY (the attack-rate profile itself differs).
The owner's test: biology should show the same sign in every country, surveillance should not; the
external anchor is the PHIRST cohort (South Africa, twice-weekly PCR irrespective of symptoms,
Cohen et al. 2021, Figure 2A), whose infection incidence by age collapses to our groups as
young/adult 1.65-1.80 and elderly/adult 0.80.

**Design.** `code/06_comp_model/run_age_experiment.R`: the 12 countries with >= 6 age-complete seasons,
four variants -- A none, B reporting offsets, C per-contact susceptibility (row scaling of the contact
matrix renormalised to spectral radius 1, so R0_s keeps its meaning), D both. Nested variants are
warm-started from their parent's optimum and polished from the parent's EKF optimum (the objective is
bimodal: a reporting-like and a susceptibility-like mode; without this D 'lost' to B by optimiser
noise). Gains are pure EKF log-likelihoods (each N(0,1) age prior adds 0.92 nats at its centre).
`run_susc_grid.R`: the profile likelihood of ONE susceptibility profile shared by all countries
(young x elderly grid) with the reporting offsets free per country -- the across-country separation.

**Findings.**
1. *Direction.* Under age-invariant reporting the contact structure alone under-predicts the young
   in 10/12 countries and the elderly in 11/12 (relative to adults; exceptions NO and IE for the
   young, HR for the elderly). The same sign almost everywhere = an age effect the model lacks; but
   the magnitudes (young 0.45-5.2x, elderly 0.95-7.2x) are far too heterogeneous for biology, so
   country-specific reporting is needed whatever else is true.
2. *Per-country likelihood.* Reporting offsets beat susceptibility in 9/12 countries (median gain
   over A: 188 vs 130 nats; |difference| > 2 nats in 11/12); susceptibility wins in EE (+56), IE (+5),
   BE (+2). Both together add a median 8 nats (-1 to +55) over the better single mechanism, and
   within a country the elderly effect wanders freely between the two -- the trade-off predicted.
3. *Attack-rate profile vs PHIRST.* With contacts alone (A, B): young/adult 0.95 (0.74-1.17),
   elderly/adult 0.29 (0.16-0.47). Per-country susceptibility (C) reproduces PHIRST's young ratio
   in the median (1.76) but ranges 0.34-6.6 across countries; elderly 0.50.
4. *Shared profile, reporting free (the decisive test).* Total over 12 countries: best at young 1.0,
   elderly 1.7 (+89 nats; elderly 2.8: +39). Any extra young susceptibility LOSES likelihood in
   almost every country (1.4: -55 to -89 total; 1.8: -210 to -290; 2.4: about -500; per country only
   PL wants 1.8 and DK is flat at 1.4). Every country's own best profile has elderly > 1 (2.8 in 8 of
   12), with gains of 1-13 nats each (EE 59). Under it the elderly/adult attack ratio is 0.5-0.9
   (BE 0.78, DK 0.86, IT 0.90, NO 0.90), close to PHIRST's 0.80; the young/adult ratio stays ~1.

**Reading.** The young excess in the observations is a LEVEL effect: once each country's reporting
level is free, the weekly wave shapes reject the young being more susceptible per contact than the
contact matrix implies (a higher young share of transmission would make their wave earlier and
sharper relative to the adults', and the data do not show it). It is reporting -- children consult
more per infection, paediatric sentinel practices -- with a magnitude that is a property of each
surveillance system. The elderly excess has a DYNAMIC component with the same sign in all 12
countries: the elderly are infected roughly 1.7-2.8x more per contact than the (vaccination-adjusted)
contact matrix implies -- immunosenescence, care homes, or contact matrices that understate the
elderly's exposure. That is the shared-biology signature the owner's test asked for, and it moves
the modelled elderly attack rate to where the reporting-independent cohort puts it.

**Unresolved.** PHIRST's young/adult infection ratio of 1.7 against the model's ~1.0 at any
plausible young susceptibility. The remaining candidate is age-specific INITIAL immunity (S0 by age,
open decision E2/11.1): adults carry more prior immunity than children, which raises the young's
attack rate with a different shape signature than per-contact susceptibility. A grid over
S0_young/S0_adult with reporting free is the next test; child-child contact weights of the synthetic
matrices are the other suspect.

**DECIDED by the owner, 2026-09-11.** Three decisions, taken on the evidence above:
1. ONE elderly susceptibility factor SHARED ACROSS COUNTRIES, with the YOUNG FIXED AT 1. Biology is
   one number for Europe, not twelve; and the grid showed a free young factor loses likelihood in 10
   of 12 countries, so fitting it would only let it re-absorb reporting level. Prior
   `log2 sigma_eld ~ N(1, 0.5)`: centre 2x, 95% band 1.0-4.0x.
2. AGE REPORTING OFFSETS ARE FITTED, one pair per COUNTRY, SHARED ACROSS that country's SEASONS --
   reporting is a property of a surveillance system, not of a season.
3. AGE-SPECIFIC INITIAL IMMUNITY IS PARKED (not rejected) to keep the model from over-complicating
   and to avoid a third age parameter trading off against the two we keep. Consequence, to be
   REPORTED rather than fitted: the model's absolute attack rate in the young is probably too low,
   and PHIRST's young/adult ratio of ~1.7 against the model's ~1.0 stands as a documented shortfall.
4. The SEASON OBSERVATION DEVIATION is ONE VALUE PER SEASON, SHARED ACROSS COUNTRIES (8 parameters,
   not 96): the observation-side twin of the shared season transmissibility. Interpretation: ILI+
   per infection in that season relative to the norm -- strain symptomaticity (viral, shared), plus
   whatever season size the mechanism cannot produce. The per-country deviations already moved
   together (FR-ES 0.97, FR-NO 0.82, DK-ES 0.78); Estonia is the known exception, so
   country-season misfit now lands in phi, S0_c or R0_s and must be watched. Per-country fits keep
   the deviation free per season as the warm start and as the diagnostic of the sharing assumption.
PHIRST remains a check, not a target. D (both mechanisms free per country) is not to be fitted --
it is unidentified within a country. Recorded in ASSUMPTIONS.md C7, F2, E5, 11.1, 11.2.

**Also this round.** Seed prior recentred to 10^-6.5 (sd 3) after 29% of fitted seeds fell below the
old prior's lower bound (late waves need small seeds); the EKF post-update clamp floored at 1e-12
instead of 0 and re-seeded the vaccinated infectious weekly (fixed in both engines, tested); the
fitter's start vector is built by one function checked for every switch combination (a missing slot
had silently shifted every later parameter); the parameters are spelled out in words in
ASSUMPTIONS.md section 12 (phi = measurement noise, does not propagate; q = process noise, does).

## 2026-09-12 compartmental model: identifiability, measured

**Method.** At each country's fitted optimum, the Hessian of the PURE negative log-likelihood (the
observed Fisher information: the curvature of the fit-quality landscape, with the priors switched off)
and of the penalised objective. Eigenvalues are the steepnesses of that landscape, eigenvectors the
directions; a near-zero eigenvalue is a parameter COMBINATION the data cannot see at all. The
difference between the two matrices is exactly what the priors hold up. Twelve countries, 86
country-seasons, current defaults (age reporting offsets on, no age susceptibility).

**Finding 1 -- an EXACT flat direction, in all 12 countries.** Every spectrum contains an eigenvalue
of magnitude 2.6e-8 to 2.3e-7 against a largest of 1.3e5 to 3.4e5, i.e. zero to machine precision.
Its eigenvector matches the reporting symmetry with cosine 1.000: the likelihood depends on the
country reporting level `c_c` and the season deviations `exp(delta_s)` only through their products, so
`c -> c*k` with `delta_s -> delta_s - log k` is numerically invisible. Nine knobs, eight products. On
Denmark 89% of the flat direction sits on the eight deviations and 11% on the country level. Today only
the weak N(0, 0.5) prior on the deviations decides where along that valley the fit stops, so the
ABSOLUTE reporting proportion is an assumption read back, not a result. FIX (required before the joint
fit): constrain the deviations to a geometric mean of one. With the deviations SHARED across countries
(decision 2026-09-11) this matters more, not less: one flat direction would then slide all twelve
country reporting levels together as a block.

**Finding 2 -- susceptibility is ENTANGLED, and the joint fit is the cure.** `S0` never appears alone:
only as `R0_s*S0` (rise rate), `c*S0` (level) and `I0_s/S0` (arrival). Conditionally -- everything else
known -- the data pin logit `S0` to sd 0.0067. Marginally, with the eight season `R0_s` free, the
likelihood-only sd is 1.04 (median over 11 countries; range 0.75-2.15), i.e. a 150-fold inflation from
the trade-off. On Denmark that is a 95% range of S0 = 0.27-0.985: effectively nothing. Treat the
`R0_s` as KNOWN and the same data give 0.825-0.839. Median sharpening 50x (range 11-160x). This is the
quantitative case for the joint fit: pooling `R0_s` over 12 countries removes the partner `S0` trades
with. CAVEAT: the conditional calculation is an UPPER BOUND, since a joint fit pins `R0_s` with finite
precision; but 8 numbers carried by 86 waves should recover much of it.

**Finding 3 -- today's S0 error bars are borrowed from the R0 prior.** The penalised fit reports logit
`S0` to sd 0.112, far tighter than the likelihood alone (1.04) or its own prior (1.0) can give. The
tight `log R0 ~ N(log 1.5, 0.05)` prior is what breaks the entanglement from outside. PREDICTION to
test: widening that prior leaves the point estimates alone (as the earlier sensitivity test found) but
inflates S0's interval substantially. If so, the width of the R0 prior is a scientific assumption that
must be argued for, not chosen for convenience.

**Per-parameter contraction (1 - sd_post/sd_prior; 0 = the posterior is just the prior).** Age
reporting offsets 0.91, reporting `c` 0.91, `S0` 0.89, season seeds 0.83, `phi` 0.80, season deviations
0.61, `R0_s` 0.60. Note the trap: `c` and the deviations look sharply contracted yet their absolute
level is not identified at all -- conditional precision and marginal precision are different questions
and only the second is publishable.

**Three code defects it exposed.**
1. `comp_model_fit.R:28` and `:122` hard-code two baseline slots (`b_RespiCompass`, `b_ERVISS`)
   regardless of which sources a country's seasons come from. ES and NO have ZERO ERVISS seasons, so
   their `b_ERVISS` enters no likelihood term and carries no prior: zero curvature from either source,
   the penalised Hessian is EXACTLY SINGULAR, and Laplace standard errors are impossible. Both needed
   the slot detected and dropped by hand. Build the slots from the sources present.
2. Where a source contributes ONE season its baseline is effectively free (data-only sd up to 2.1e3 on
   the log scale). `b` is the only unpenalised parameter, so nothing catches it. CZ, IT, NL each have a
   single ERVISS season. A weak prior on `b` would close this.
3. PL (1 direction) and NL (2) have NEGATIVE likelihood curvature at the optimum. Legitimate rather
   than a convergence failure -- the optimum minimises the PENALISED objective, so the prior may supply
   the missing curvature -- but the prior is doing structural work there, and PL's `S0` inverse was
   numerically degenerate as a result.

**Visual summary.** Parameter scope map plus these results:
https://claude.ai/code/artifact/f6fddc68-046e-4d18-9d75-ce896edaf16f

## 2026-09-12 the working model: simpler, joint, and recovery-tested

**Decision (owner).** This is the WORKING MODEL. More iterations are expected -- "this is how science
works" -- but it is the base everything else builds on. Spain is included (the season cut moved from 6
to 5, giving 12 countries, 86 country-seasons, 184 parameters). Uncertainty intervals are required
output. The Kalman filter and the driver/subtype validation both WAIT.

**What it is.** `code/07_joint_model/`, with the one-paragraph abstract and the full cut list in
MODEL.md. A joint age- and vaccination-structured SIR over 12 countries and 8 seasons: transmissibility
varies between seasons and is shared across countries; susceptibility varies between countries and is
shared across their seasons; one global elderly susceptibility; reporting varies by country and age but
not by season except for one shared season visibility constrained to average one; a free seed per
country-season for arrival time; negative-binomial counts with a country dispersion.

**What was dropped, and why it was defensible.** The Kalman filter, and with it the process noise and
the initial covariance: measured on the pilot, the filter only behaved with both pinned, and with them
pinned its optimum equalled the deterministic one in every parameter. Dropping it also freed the
observation model from the Gaussian a Kalman update requires, which mattered because 29% of observed
cells are exactly zero and the pilot's fitted dispersion put 15% of its predictive mass below zero.

**Architecture.** The entire log-posterior is one C++ call, so R does no per-evaluation work (the pilot
spent 42% of each objective in R glue and 28% building trajectories the optimiser discarded). No filter
means no Jacobian, which is most of why the core is short. Fitting exploits separability: only 16 of 184
parameters are shared across countries, so block coordinate descent optimises all 12 country blocks in
parallel, then the shared block, then polishes jointly. 116 s for the full fit; the C++ is ~1750x the
base-R reference it is tested against.

**Two optimiser findings that cost real time.**
1. Each country's local block has a SECOND, WRONG OPTIMUM: let the dispersion collapse and the negative
   binomial becomes so diffuse that every curve fits, so nothing forces a wave and the country flat-lines.
   The first joint fit lost the Netherlands that way. It is a LOCAL OPTIMUM, not a prior problem --
   tightening the dispersion prior moved the failure to Poland instead. Fixed by multi-starting each
   local block from three points; worth 406 nats of likelihood and took the unidentified directions from
   9 to 0. A residual 1-2% per-country failure rate remains and must be checked for on every fit.
2. Block coordinate descent leaves a slow tail when the shared and local blocks are correlated (about 2
   nats per sweep at the cap). The joint polish recovers it (24-31 nats), so keep both stages.

**Identifiability, measured.** Every parameter family contracts between 0.69 and 0.94 against its prior,
so nothing in this model is a restatement of an assumption. Contrast the pilot, where susceptibility
looked well determined only because it was borrowing the tight transmissibility prior.

**The honest limitation.** Every country needs about 2.6x more observation noise than its own
week-to-week scatter can explain (`jm_adequacy`, figure 05). That excess is the deterministic mean
failing to follow the wave, written off as measurement error. It is the trigger condition for restoring
the filter, and the filter should be judged on whether it CLOSES THAT GAP rather than on whether it
moves the estimates.

**What it learns, with intervals.** Season visibility spans 0.57-1.80 with non-overlapping intervals
between the extreme seasons, while transmissibility spans 1.51-1.71 with intervals of about +/-0.09 that
mostly overlap. So between-season differences in observed burden are mostly about how VISIBLE a season
was, not how TRANSMISSIBLE it was. That is the pilot's Danish conclusion, now carried by 86 waves.
Caveat to carry: season visibility is partly a residual absorber, so read it with the noise budget.

**Recovery.** See MODEL.md for the numbers. Headline: the publishable quantities recover (rank 0.95-0.97,
coverage 96% median), the intervals are conditional on the optimiser finding the right basin (coverage
falls to 25-79% when truth is drawn from the priors), and the learning layer's two-step procedure is
UNBIASED, so driver slopes can be read at face value.

**The recovery harness is built to carry the learning layer.** `jm_truth_with_driver` constructs a world
where a covariate really moves the season parameters by a known amount, and `jm_driver_recovery` runs the
fit-then-regress step and checks the slope comes back. Any driver association found in real data can
therefore be told apart from one manufactured by the procedure, before it is ever claimed.
