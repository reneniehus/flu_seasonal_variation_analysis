# The joint EU/EEA influenza model, version 1

One page. What it assumes, why it is simpler than the compartmental pilot, and what was cut.

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

| Parameter | Varies | Count (12 countries, 8 seasons) | Meaning |
|---|---|---|---|
| `R0_s` | between seasons, same across countries | 8 | transmissibility of that season's virus at full susceptibility |
| `delta_s` | between seasons, same across countries, averages zero | 8 (7 free) | that season's positive consultations per infection, relative to normal |
| `sigma_eld` | nothing: one global value | 1 | elderly susceptibility per contact, relative to adults |
| `S0_c` | between countries, same across seasons and ages | 12 | susceptible fraction at season start |
| `c_c` | between countries, same across seasons | 12 | positive consultations per infection |
| `off_c,young`, `off_c,eld` | between countries and age groups | 24 | age offsets on reporting, adults the reference |
| `phi_c` | between countries | 12 | negative-binomial dispersion of the weekly counts |
| `b_c,src` | between countries and data sources present | 22 | off-season floor, per data source |
| `I0_c,s` | freely between country-seasons | 86 | seed size, i.e. arrival time of that wave |

Total about 184, of which 16 are shared across countries and 168 are local to one country. The
likelihood is therefore separable given the shared block, which is what the fitting strategy exploits.

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
fixed, all twelve in parallel, then the shared block is optimised with the locals fixed, repeated to
convergence and finished with a joint polish. That exploits the separability above: a country's own
likelihood is a twelfth of the joint one, so a sweep costs about 33 full-likelihood evaluations against
185 for one finite-difference gradient of the flat problem.

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

**Residual defect.** In 1 of 8 replicates, 1 of 12 countries fell into the flat-line optimum
(dispersion collapses, no wave fitted): roughly a 1-2% failure rate per country-fit. The multi-start
of the local blocks reduced this from 1-in-11 but did not remove it. It is DETECTABLE rather than
silent -- a flat-lined country is obvious in figure 03 and in the per-country correlations -- so
check for it on every real fit.
