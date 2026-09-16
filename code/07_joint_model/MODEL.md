# The joint EU/EEA influenza model, version 1

One page. What it assumes, why it is simpler than the compartmental pilot, and what was cut.

> ## ⚠ ONE OPEN QUESTION ABOUT THE SOURCE DATA — it affects what the model is fitted to
>
> **Does an ERVISS week with `tests > 0` and no `detections` row mean zero detections, or an
> unpublished count?** 3151 weeks state `detections = 0` explicitly; 709 omit the row with tests
> recorded, and no published positivity value rescues any of them.
>
> We re-derive positivity as detections/tests, so the first becomes an observed zero and the second
> becomes missing — and a missing week is dropped, which for a trailing run shortens the season.
>
> **Pending an answer we exclude the country-seasons containing weeks that CANNOT plausibly have been
> zero** — judged per week from that week's test count against the local positivity of its published
> neighbours, not from a count of affected weeks. In the fitted design that is **CZ 2024/2025 alone**
> (14 of its 17 affected weeks are not plausibly zero, 13 of them with no published neighbour at all,
> because CZ's detections feed went dark on 2025-03-26). **PL 2024/2025 is kept**: its single affected
> week sits among neighbours with 0 detections over 114 tests, so it is almost certainly a real zero.
>
> The design is therefore **12 countries, 8 seasons, 85 country-seasons, 182 parameters**. All 12
> countries and all 8 seasons survive. CZ also loses its ERVISS baseline slot, because 2024/2025 was
> its only ERVISS season — the right outcome, since that slot was otherwise fitted with no off-season
> behind it.
>
> **The full question, the evidence, and what changes under either answer are in
> `documentation/to_confirm_with_surveillance.md` (question 1)** — along with five other things we
> have inferred about the source data rather than confirmed. `erviss_encoding_ambiguous()` in
> `stitch_iliplus.R` is the implementation; `jm_build_data(exclude_ambiguous_positivity = FALSE)` fits
> everything anyway.

## Abstract

We fit weekly influenza-positive ILI consultations from twelve EU/EEA countries over eight seasons
with a single age- and vaccination-structured SIR, estimated jointly. Three age groups (0-14, 15-64,
65+) mix by each country's contact matrix, rescaled so its dominant eigenvalue is one, which makes a
season's basic reproduction number the realised R0 given that mixing rather than an artefact of the
contact data's absolute scale. Each country-season is an independent epidemic, seeded on 1 August and
integrated on a daily grid, with a single vaccination pulse into the 65+ group on 1 October. What
varies where is the model's central design and is chosen so that each quantity is identified by a
different feature of the data. Transmissibility `R0_s` is assumed to **vary between seasons and to be
the same across countries**, because the dominant subtype is a Europe-wide property of the virus;
initial susceptibility `S0_c` is assumed to **vary between countries and to be the same across that
country's seasons**; and the elderly's relative susceptibility per contact is **assumed the same
everywhere**, one number, because it is biology rather than surveillance. On the observation side the
reporting proportion, the expected positive consultations per infection, is assumed to **vary between
countries** and, through two offsets, **between age groups within a country**, but to be **the same
across that country's seasons** apart from one Europe-wide deviation per season that is **constrained
to average zero**, so the country level and the season deviation are separately identified rather than
sliding against each other. The seed size is left to **vary freely between country-seasons**: it is
the model's only handle on when a wave arrives, and it cannot be shared without forcing every season
to arrive at the same time. Observed counts are negative-binomial around the deterministic mean with a
country-specific dispersion, which lets the 29% of weeks that read exactly zero carry their proper
weight instead of being pushed through a Gaussian that would place a sixth of its mass below zero. The
infectious period, the three vaccine effects and the vaccination timing are fixed from external
estimates; no parameter is allowed to vary in more than one direction at a time.

## What was cut relative to the compartmental pilot, and why

| Cut | Why |
|---|---|
| **The Kalman filter, and with it the process noise `q` and the initial-state covariance `P0`** | Measured: the filter only behaved when `q` was fixed small and `P0` set to zero, and with those the EKF optimum agreed with the deterministic one in every parameter (S0 0.830 vs 0.834, `R0_s` within 0.02). It was costing 7x the runtime and two unfittable knobs to reproduce the answer the deterministic fit already gave. Retained as an option to add back once the innovation diagnostics exist. |
| **The Gaussian observation with count-like variance** | Only needed because a Kalman update requires a Gaussian. Dropping the filter frees the natural choice. At the fitted dispersion the Gaussian put 15% of its predictive mass below zero at any count level, and 29% of observed cells are exactly zero. |
| **The two-stage fit, the engine switch, the age-susceptibility switch, the per-country `R0_s` diagnostic stage** | Machinery that existed to manage the filter or to answer questions now answered (decisions.md, 2026-09-11). |
| **Per-country-season season deviation** | Now one shared value per season (owner decision), which is both the scientific hypothesis and 88 fewer parameters. |
| **Age-specific initial immunity, age-specific susceptibility of the young** | Parked by decision. The young's absolute attack rate is therefore probably too low; recorded as a limitation, not fitted. |

## What was fixed relative to the pilot

- **The exact flat direction is gone.** The season deviations are constrained to average zero on the
  log scale. Without it the likelihood is invariant to raising the country reporting level and lowering
  every deviation by the same factor, a direction measured at eigenvalue 1e-7 against 1e5 in all twelve
  countries, and with the deviations now shared that one direction would slide all twelve country
  levels together.
- **Susceptibility is un-entangled by construction.** `R0_s` is shared, so it is pinned by 86 waves and
  can no longer absorb a country's `S0`. Measured gain from the pilot: a 50-fold sharpening in the
  likelihood-only precision of `logit S0`.
- **Baseline slots are built from the sources a country actually has**, and carry a weak prior. In the
  pilot two slots were created unconditionally; Spain and Norway have no ERVISS seasons, so their
  second slot had zero curvature from data or prior and the Hessian was exactly singular.
- **The `R0_s` prior is widened** from sd 0.05 to 0.15. In the pilot the tight prior was what broke the
  susceptibility entanglement, so the reported S0 intervals were largely inherited from it. Sharing
  `R0_s` does that job with data instead.

## Parameters

| Parameter | Varies | Count (12 countries, 8 seasons, 85 country-seasons) | Meaning |
|---|---|---|---|
| `R0_s` | between seasons, same across countries | 8 | transmissibility of that season's virus at full susceptibility |
| `delta_s` | between seasons, same across countries, averages zero | 8 (7 free) | that season's positive consultations per infection, relative to normal |
| `sigma_eld` | nothing: one global value | 1 | elderly susceptibility per contact, relative to adults |
| `S0_c` | between countries, same across seasons and ages | 12 | susceptible fraction at season start |
| `c_c` | between countries, same across seasons | 12 | positive consultations per infection |
| `off_c,young`, `off_c,eld` | between countries and age groups | 24 | age offsets on reporting, adults the reference |
| `phi_c` | between countries | 12 | negative-binomial dispersion of the weekly counts |
| `b_c,src` | between countries and data sources present | 21 | off-season floor, per data source (CZ has only RespiCompass once its single ERVISS season is excluded) |
| `I0_c,s` | freely between country-seasons | 85 | seed size, i.e. arrival time of that wave |

Total **182**, of which 16 are shared across countries and 166 are local to one country. The
likelihood is therefore separable given the shared block, which is what the fitting strategy exploits.
(Without the provisional positivity-encoding exclusion it is 86 country-seasons and 184 parameters —
`jm_build_data(exclude_ambiguous_positivity = FALSE)`; both designs are pinned by the test suite.)

## What the data layer does and does not support (audited 2026-09-16)

The model code had been audited; the data layer under it had not. It reconstructs correctly — all 86
country-season observation matrices rebuild **exactly** from the raw streams by independent code (max
count difference 0), the age bands and populations are mutually consistent to the person, missing cells
are genuinely missing, and the season boundaries are the documented 1 August. Five things nonetheless
qualify what can be read off the output (a sixth, the positivity-encoding
exclusion, is in the box at the top of this file and in
`documentation/to_confirm_with_surveillance.md`):

- **Season support is uneven, and the thinnest season carries the highest `R0`.** `2025/2026` rests on
  **7 of the 12 countries** (DK, EE, FR, BE, IE, PL, HR), on grids of 34–41 weeks against 52–53
  elsewhere, and 774 observed cells against 1497 for 2023/2024 — and its fitted `R0` of 1.71 is the
  highest of the eight. `jm_summary_season` now reports `n_country`, `obs_cells` and the grid range
  alongside `R0`, so no season-level number is read without its sample size.
- **Norway has no contact matrix of its own** and is fitted on the EU average collapsed with Norwegian
  populations. That matters because the age reporting offsets are the parameters most sensitive to the
  mixing pattern: re-optimising Norway's block under neighbouring countries' matrices moves its elderly
  offset by up to 66% for a likelihood spread of 0.7 nats. `d$contact_source` and
  `summary_country.csv` now name each country's matrix.
- **The attack rate is now integrated to a fixed 53-week horizon** for every country-season. It used to
  stop wherever that country's surveillance series stopped, so it was a window quantity reported as a
  season quantity: five cells moved by more than 5% and IT 2015/2016 by 15% (window 38 weeks) for no
  epidemiological reason. The likelihood is unchanged by this — verified to the last bit.
- **A country-season's baseline can be fitted partly on the other source's weeks.** The panel stitches
  per week, so 8 of the 86 (BE, CZ, DK, EE, HR, IE, NL, PL 2023/2024) take 3–14 weeks from ERVISS while
  labelled RespiCompass, and for 7 of the 8 those weeks sit in the late off-season tail where the
  baseline is essentially the whole model mean. The stitch documents the mixing; its consequence for
  `b_c,src` is this: read `b` as a coverage-era nuisance parameter, not as a property of a surveillance
  system.
- **The twelve countries are a deliberate selection, not what the data admit.** At the design's own
  inclusion rule, Iceland (7 age-complete seasons), Malta (6) and Austria (5) would also qualify and
  were never offered to the model. Since sharing `R0_s` across the countries in the design is what
  identifies `S0_c`, which countries are in it is a substantive choice.

## Fixed, not fitted

`gamma = 1/3.6` per day; `ve_inf = 0.25`, `ve_ili_cond_inf = 0.20`, `ve_spread = 0.20`; vaccination
pulse on season day 62 into 65+ only, at the reported national coverage; the contact matrix, rescaled
to spectral radius one; susceptibility shared across age groups; no waning, no ageing, no importation
after the seed, no latent period, single strain.

## Implementation

`joint_model.cpp` holds the whole log-posterior: the daily SIR for every country-season, the
negative-binomial likelihood, and the priors. One call from R returns one number, so no per-evaluation
R overhead. `joint_model.R` assembles the data, packs and unpacks the parameter vector and drives the
fit by **block coordinate descent**: each country's local block is optimised with the shared block
fixed, all twelve in parallel, then the shared block is optimised with the locals fixed, repeated and
finished with a joint polish. That exploits the separability above: a country's own likelihood is a
twelfth of the joint one, so a sweep costs about 33 full-likelihood evaluations against 185 for one
finite-difference gradient of the flat problem. 100-170 s for the full 12-country fit depending on
machine load; the curvature intervals and the two Hessians behind the identifiability report cost
another 7 minutes or so, which is why a recovery replicate that does not need intervals skips them.

**How to read the convergence report.** The sweep loop is EXPECTED to run out of sweeps: its gain
decays geometrically (48.5, 10.2, ..., 2.5, 2.3 nats) and never reaches a 0.05-nat tolerance in 15,
extrapolating to a tail of about 24 nats. The joint polish is what clears that tail, and it gained 31
nats on the real fit -- slightly more than the extrapolation, the difference being cross-block
curvature the sweeps cannot see by construction. So the sweep loop's status is reported separately
(`sweeps_hit_tol`, `sweep_tail_nats`) and `converged` means what a reader expects: all twelve local
blocks, the shared block and the joint polish all returned success, and no country is stuck flat. It
previously meant only "the sweep loop did not run out of sweeps", so it read FALSE on a fit where all
three stages had succeeded -- a flag that is always FALSE is a flag nobody reads.

Two things guard the optimiser, both of them there because the failure they prevent actually happened.

**Multi-start on every local block.** Each country's block has a second, WRONG optimum: let the
dispersion collapse and the negative binomial becomes so diffuse that every curve fits about equally
well, so nothing pushes the model to place a wave and the country flat-lines at its baseline. The
first joint fit lost the Netherlands that way. It is an optimiser failure, not a fact about the data --
a harder search found a solution 406 nats better. Each block is therefore started from three points
(the warm start, the same with the dispersion pinned at the value that country's own week-to-week
scatter implies, and that plus an earlier arrival) and the best is kept. The block is also seeded with
the value it came in with, so it can never be written back worse.

**The flat-line protector** (`jm_flat_check`, `jm_unflatten`) runs inside `jm_fit`, before and after
the polish, so the returned fit is always checked. Detection is MECHANISTIC rather than
goodness-of-fit based, because a country can fit badly for honest reasons but cannot have an epidemic
that never happened: the trigger is the attack rate collapsing (threshold 0.03 against a measured
minimum of 0.30 across the twelve countries) or the fitted peak being almost all baseline (0.15
against a measured 0.99). Per-season columns are reported too, and a PARTIAL flat line triggers when
one season is flat while at least two others are healthy -- without that, three of the Netherlands'
six seasons could flat-line invisibly at a cost of 3084 nats. A rescue is adopted only when it
improves that country's objective; if the flat solution still wins after the whole escape ladder it is
reported as UNRESOLVED rather than papered over, because that means the series carries no identifiable
wave and forcing one would hide a finding about the data. Residual risk, measured on the recovery
replicates: roughly 1-2% per country-fit, and detectable rather than silent.

## Does it recover a known truth? (measured 2026-09-12)

`run_joint_recovery.R`, 18 full refits on the real 12-country design. Data simulated from the model
itself at a known parameter set, then the whole pipeline refitted from scratch.

**1. Truth at the fitted optimum, 8 replicates.** The publishable quantities come back:

| quantity | rank recovery (Spearman) | 95% interval coverage |
|---|---|---|
| `R0_s` by season | 0.95 | 98% |
| season visibility `delta_s` | 0.96 | 89% |
| `S0_c` ranking across countries | 0.97 | 99% |
| `sigma_eld` | -- | 100% |

Median coverage across families 96% against a nominal 95%, so on data resembling ours the intervals
mean what they say.

**2. Truth drawn from the priors, 4 replicates -- the warning.** Rank recovery stays high (`R0_s`
Spearman 0.99) but COVERAGE COLLAPSES to 25-79% and root-mean-square error rises by an order of
magnitude. The fit occasionally lands in a different optimum on harder configurations, and the
interval is then the wrong width in the wrong place. Consequence for reporting: the intervals are
calibrated for data like ours and are CONDITIONAL ON THE OPTIMISER HAVING FOUND THE RIGHT BASIN; they
are not a general guarantee.

**3. Driver recovery, 6 replicates -- the learning layer's own validity test.** Truth built so a
covariate really moves the season parameters (`jm_truth_with_driver`), then the two-step
fit-then-regress procedure run on the FITTED season parameters (`jm_driver_recovery`):

| effect | truth | recovered (mean) | sd across replicates |
|---|---|---|---|
| covariate on `log R0_s` | 0.060 | 0.059 | 0.009 |
| covariate on `log delta_s` | 0.350 | 0.334 | 0.074 |

Essentially UNBIASED, both within one sd of truth. So the learning layer may regress fitted season
parameters on drivers and read the slopes at face value; no attenuation correction is needed. (An
earlier single replicate on a cut-down 2-country design suggested ~20% attenuation. That was an
artefact of n = 1 on the wrong design and does not hold.)

**4. Can the recovery test fail at all? (measured 2026-09-14) -- the sharpest limitation.** Studies
1-3 all simulate from the model's own parameter space, so on their own they ask whether the OPTIMISER
works and would certify the model however wrong its assumptions are about the world. So the harness
also simulates from truths the model CANNOT represent (`jm_simulate_violation`,
`jm_misspecification_check`), one replicate per arm on the real 12-country design:

| arm | `R0_s` rank | visibility rank | **`S0_c` rank** | noise excess |
|---|---|---|---|---|
| control, truth representable | 1.00 | 0.96 | **0.97** | 1.25x |
| each country's own `R0` (sd log 0.10, i.e. spanning 0.84-1.20) | 0.95 | 0.96 | **0.09** | 1.26x |
| a second wave the SIR cannot make (+60% bump late in every season) | 0.93 | 0.96 | **0.99** | 1.24x |

Three things follow, and all three matter more than any number in studies 1-3.

- **The test is not vacuous.** It fails, decisively, on the project's target quantity.
- **The season-level conclusions are robust; the country ranking is conditional.** Transmissibility and
  visibility recover under both violations. But `S0_c` is identifiable ONLY BECAUSE `R0_s` is shared
  across countries, and that is a double edge: if countries truly differ in transmissibility, the
  difference has nowhere to go but into `S0_c`, and the susceptibility ranking degrades to noise.
  A 10% between-country spread in `R0` is entirely plausible (school calendars, climate, housing), so
  this is not an extreme stress test. **Report the `S0_c` ranking as conditional on shared
  transmissibility, and say so.**
- **The noise budget cannot see either violation** (1.24-1.26x against the control's 1.25x), so it is
  not the diagnostic for this. Note also what that implies about the real-data 2.59x: a 60% spurious
  late wave in every season costs essentially nothing in the noise budget, so the real gap is not a
  missed secondary wave but something pervasive across the whole season -- which is what process noise
  describes and what the filter should be judged against.

The diagnostic that WOULD settle it is a per-country `R0` multiplier fitted as an alternative model and
compared by likelihood. That is a model comparison, which is the learning layer's first job, and it is
now the top item on that list.

**Residual defect.** In 1 of 8 replicates, 1 of 12 countries fell into the flat-line optimum
(dispersion collapses, no wave fitted): roughly a 1-2% failure rate per country-fit. The multi-start
of the local blocks reduced this from 1-in-11 but did not remove it. It is DETECTABLE rather than
silent -- a flat-lined country is obvious in figure 03 and in the per-country correlations -- so
check for it on every real fit.

## Are the headline claims robust? (stress-tested 2026-09-14)

Two results carry the model's weight, so both were attacked directly rather than reported as found.

**1. "Seasons differ more in visibility than in transmissibility" is the data, not the priors.** The
worry is real: the priors are asymmetric by construction (`log R0 ~ N(log 1.5, 0.15)` against
`delta ~ N(0, 0.5)`), so a tight prior on one and a loose prior on the other could manufacture
exactly the asymmetry we report. Refitting under four prior settings says it does not:

| prior setting | transmissibility sd(log) | season visibility sd(log) |
|---|---|---|
| as fitted (R0 0.15, visibility 0.50) | 0.036 | 0.377 |
| R0 prior widened 4x to 0.60 | 0.037 | 0.377 |
| visibility prior tightened to 0.15 | 0.038 | 0.322 |
| **fully symmetric, both 0.50** | **0.037** | **0.377** |

Visibility varies about ten times as much as transmissibility under every setting, including the
symmetric one, and widening the `R0_s` prior fourfold changes its spread by 0.001. What the prior
*does* control is the LEVEL: widening it moved the fitted `R0_s` range from 1.51-1.71 to 1.54-1.76.
So report the spread as a finding and the level as prior-informed.

**2. The 2.6x noise gap is not an artefact of how the scatter is measured.** The comparison needs a
choice of smoothing window and of which weeks count as epidemic, and either could be doing the work.
Across 18 specifications (3-, 5- and 7-week moving average; epidemic threshold 5, 20 and 50 per
100 000; with and without the small-sample correction) the excess runs **2.13x to 3.65x**. The
reported 2.59x sits in the lower half of that range, so the gap is robust and the headline number is
the conservative end of it, not the flattering one.
