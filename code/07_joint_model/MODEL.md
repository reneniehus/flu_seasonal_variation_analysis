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
> The design is therefore **12 countries, 8 seasons, 85 country-seasons, 181 parameters**. All 12
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

We fit weekly influenza-positive ILI consultations from **twelve EU/EEA countries over eight seasons —
85 country-seasons, 8 685 observed age-week cells — with a single age- and vaccination-structured SIR,
estimated jointly in 181 parameters**. **R0 is fixed at 1.5 from the literature, not fitted.** Three age groups (0-14, 15-64, 65+) mix by each country's
contact matrix, rescaled so its dominant eigenvalue is one, which makes a season's basic reproduction
number the realised R0 given that mixing rather than an artefact of the contact data's absolute scale.
The elderly susceptibility re-weights that matrix and the result is **rescaled again**, so `sigma_eld`
redistributes *who* gets infected without changing *how transmissible* the season is: `R0_s` keeps its
meaning whatever `sigma_eld` does, and `sigma_eld` is therefore identified by the age composition of
cases rather than by the size of the wave.
Each country-season is an independent epidemic, seeded on 1 August and integrated on a daily grid to a
fixed 53-week horizon, with a single vaccination pulse into the 65+ group on 1 October.

**What varies where is the model's central design**, and it is forced by the data rather than chosen
for convenience. Observed ILI+ levels differ more than a hundredfold between countries and each country
keeps its rank across seasons, so the level has to be a **country** property — a reporting proportion —
not an epidemiological one. Seasons rise and fall **together** across Europe once that country level is
divided out, which is what licenses one Europe-wide season effect on each quantity. Transmissibility
and susceptibility enter a wave's rise rate as a **product** and are indistinguishable from one wave,
so one of them is pinned: `R0` is fixed and **susceptibility is the single sensor of how easily a
season spread**. It is decomposed additively on the logit scale into a **country level** `S0_c`
(varying between countries, the same across that country's seasons) and a **season effect** `x_s`
(varying between seasons, the same across countries, constrained to average zero). The elderly's
relative susceptibility per contact is **assumed the same everywhere**, one number, because it is
biology rather than surveillance.

The reporting proportion has exactly the same structure on the log scale: a **country level** `c_c`,
a **season effect** `delta_s` constrained to average zero, and two age offsets **varying between
countries** because how readily children and the elderly are seen is surveillance, not biology. The
average-zero constraints are not cosmetic: reporting level and season visibility enter the observation
as a **product** and are the same number to the data, so the constraint is what separates them. The
question the design must answer from shape alone — *a bigger season: more susceptible, or more
visible?* — it can: a susceptibility effect makes the wave rise faster and peak earlier, a visibility
effect only scales it, and with `R0` fixed the rise rate is `S0`'s alone to explain. **That is why
this is a joint fit rather than 85 separate ones.** The seed size is left to **vary freely between
country-seasons**: it is the model's only handle on when a wave arrives, and the observed peak spans 14
weeks across the panel, so it cannot be shared.

Observed counts are negative-binomial around the deterministic mean with a country-specific dispersion,
which lets the **31% of observations that read exactly zero** carry their proper weight instead of being
pushed through a Gaussian that would place a sixth of its mass below zero. The infectious period, the
three vaccine effects and the vaccination timing are fixed from external estimates; no parameter is
allowed to vary in more than one direction at a time.

**What it delivers, and what it does not.** The fit converges in about 95 s; every parameter family
contracts against its prior by 0.69-0.95, and the season effect on susceptibility by 0.93, so nothing
reported is a restated assumption. The question the design must answer from shape alone -- *a bigger
season: more susceptible, or more visible?* -- it does answer: the posterior correlation between the
two season effects is only 0.27 (max 0.31 over seasons). Two limits are measured rather than asserted:
the deterministic mean needs about 2.6x more observation noise than the data's own week-to-week
scatter can explain, so some of what the model calls measurement error is really misfit; and the
country susceptibility **ranking is conditional on transmissibility genuinely being equal across
countries** -- on simulated data where countries' true R0 differed by 10% that ranking fell from 0.93
to 0.01, because with R0 pinned a real transmissibility difference has nowhere to go but `S0_c`.
Season-level conclusions survive that test; the country ranking does not. And a **model comparison**
says where the season variation prefers to live: the previous model, with a *fitted* per-season R0
and a country-only S0, fits 11.9 nats better at one more parameter (AIC +21.8 against this model).
The data lean towards seasons differing in how *fast* a wave rises rather than in how *many* are left
to infect -- a finding to carry into the learning layer, not a reason to undo the decision.

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
- **Susceptibility is un-entangled by construction.** `R0` is fixed, so nothing trades against `S0` in
  the rise rate: each wave's early growth identifies its `S0_{c,s}` outright. (The earlier joint model
  achieved this by sharing a fitted `R0_s` across countries — a 50-fold sharpening over the pilot;
  fixing `R0` keeps that identification with one parameter fewer and the founding caveat made explicit.)
- **Baseline slots are built from the sources a country actually has**, and carry a weak prior. In the
  pilot two slots were created unconditionally; Spain and Norway have no ERVISS seasons, so their
  second slot had zero curvature from data or prior and the Hessian was exactly singular.
- **No transmissibility prior at all.** In the pilot a tight `R0_s` prior was what broke the
  susceptibility entanglement, so the reported `S0` intervals were largely inherited from it. With `R0`
  pinned there is nothing to inherit: `S0`'s intervals come from the rise rates.

## Parameters

| Parameter | Varies | Count (12 countries, 8 seasons, 85 country-seasons) | Meaning |
|---|---|---|---|
| `S0_c` | between countries; the country's level at the average season | 12 | susceptible fraction on 1 August, logit scale |
| `x_s` | between seasons, same across countries, averages zero | 8 (7 free) | season effect on susceptibility, added to every country's logit `S0_c` |
| `c_c` | between countries; the country's level at the average season | 12 | positive consultations per infection, log scale |
| `delta_s` | between seasons, same across countries, averages zero | 8 (7 free) | season effect on visibility, added to every country's log `c_c` |
| `sigma_eld` | nothing: one global value | 1 | elderly susceptibility per contact, relative to adults |
| `off_c,young`, `off_c,eld` | between countries and age groups | 24 | age offsets on reporting, adults the reference |
| `phi_c` | between countries | 12 | negative-binomial dispersion of the weekly counts |
| `b_c,src` | between countries and data sources present | 21 | off-season floor, per data source (CZ has only RespiCompass once its single ERVISS season is excluded) |
| `I0_c,s` | freely between country-seasons | 85 | seed size, i.e. arrival time of that wave |

So the two quantities the project reads are each a **two-way additive decomposition**: on the logit
scale `logit S0_{c,s} = S0_c + x_s`, and on the log scale `log c_{c,s} = c_c + delta_s`. The country
levels carry the mean (their priors are centred on the plausible level, not on zero); the season
effects are centred on zero and constrained to average zero, which is what makes the level and the
effects separately identified rather than sliding against each other.

Total **181**, of which 15 are shared across countries and 166 are local to one country. The
likelihood is therefore separable given the shared block, which is what the fitting strategy exploits.
(Without the provisional positivity-encoding exclusion it is 86 country-seasons and 183 parameters —
`jm_build_data(exclude_ambiguous_positivity = FALSE)`; both designs are pinned by the test suite.)

## R0 in a country: what the contact matrix, the age susceptibility and the pyramid do, and do not do

**R0 is fixed at 1.5 in every country and every season** (owner decision, 2026-09-25). This section is
the honest account of what that means once age mixing enters, because it is easy to believe the
model has a country-specific R0 when it does not.

**How the dynamics are built.** Each country has its own contact matrix `C_c` (Prem et al., collapsed
to three groups with that country's population weights, so the pyramid already shapes it). It is
rescaled to spectral radius one. The elderly susceptibility `sigma_eld` then multiplies the elderly
row, and the result is **rescaled to spectral radius one again**. The force of infection uses
`beta = R0 x gamma` on that matrix. Consequence: the realised reproduction number at full
susceptibility is **exactly 1.5 in every country**, whatever its matrix, its pyramid, or `sigma_eld`.

**What the matrix, `sigma_eld` and the pyramid therefore control.** *Who* gets infected, not *how
many*. They set the age distribution of infections and the age-specific attack rates — the elderly
attack rate is 21% against 33% for adults precisely because of this structure — but they cannot move
the overall speed of the epidemic. That is delivered entirely by `S0_{c,s} x R0`, and with `R0`
pinned, by `S0_{c,s}` alone.

**The consequence the owner asked about, stated plainly.** In reality a country with a larger share of
elderly, who are ~2x as susceptible per contact, or with denser mixing, *has* a higher effective
transmissibility. The model erases that difference by renormalising, so wherever it is real it has
**nowhere to go but the country's `S0_c`**. `S0_c` is therefore a composite: genuine susceptibility
plus whatever transmissibility differences the demography and mixing would have produced. The same
holds across seasons for `x_s`: a genuinely more transmissible strain (H3N2 over H1N1) lands in `x_s`
as "more susceptible". This is the project's founding caveat (`decisions.md`, "Fix R0 from the
literature") and it is now the design rather than an aside. **Read `S0_c` and `x_s` as "how easily
influenza spread here / this season", not as a pure immunity measure.**

**On the visibility side the age structure IS country-specific, correctly.** Age-specific reporting
(`off_c,young`, `off_c,eld`) varies by country, and the pyramid enters the expected count through the
group populations `N_a`, so a country with more elderly sees more elderly consultations. What the
model does **not** have is age-specific *season* visibility: `delta_s` scales every age group by the
same factor, so an H3N2 season that is disproportionately visible in the elderly is represented only
in aggregate. A natural extension, not a defect.

**The alternative, for a later model comparison.** Skip the second renormalisation. Then the realised
R0 in country `c` is `1.5 x rho(D_sigma C_c)`, which varies by country through *known* demography and
mixing — a mechanistically-determined transmissibility — leaving `S0_c` a purer susceptibility
residual. It is a one-line change in the C++ and it changes what `sigma_eld` means (raising it would
then raise overall transmissibility too). Worth fitting side by side once the learning layer starts;
not done now, so that one thing changes at a time.

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

## Identifiability: the theory, then what the data say (fixed-R0 model, 2026-09-25)

**Theory.** With `R0` pinned, each wave's early growth rate `gamma (R0 S0 - 1)` identifies its
`S0_{c,s}` outright -- nothing trades against it in the rise. Its final size then follows from
`S0_{c,s}` (the SIR final-size relation), so the peak height identifies `c_{c,s}` given the shape; the
seed identifies timing. The two-way decompositions `S0_c + x_s` and `c_c + delta_s` are ordinary
additive designs on 85 cells with 12 + 7 free effects each -- 66 residual degrees of freedom -- and
the average-zero constraints fix the level. The one question the design has to answer from shape
alone is whether a bigger season is more susceptible (`x_s` up: faster rise, earlier peak, bigger
final size) or more visible (`delta_s` up: the same curve scaled). Figure 03 shows the two at equal
peak height; they differ by about five weeks in peak timing. The one exact tie that remains is
`c_c x delta_s`, broken by the constraint, not by the data.

**What the data say** (`output/joint_model/joint_fit.rds`, penalised Hessian):

| claim | the number |
|---|---|
| the data decide, not the priors | contraction 0.95 (`S0_c`), 0.93 (`x_s`), 0.94 (`c_c`), 0.90 (`delta_s`); minimum over all families 0.69 (`sigma_eld`) |
| susceptible-vs-visible is separable within a season | posterior corr(`x_s`, `delta_s`), same season: median 0.27, max 0.31 |
| country level vs reporting level is separable | posterior corr(`logit S0_c`, `log c_c`): median magnitude 0.26 |
| `S0` is read off the rise rate | Spearman(observed early growth rate, fitted `S0_{c,s}`) = 0.54 over 83 waves; the empirical rate is a crude 6-week slope, so this confirms the direction rather than the precision |
| no direction the data cannot see | 0 near-flat eigenvalues of the likelihood-only Hessian; penalised Hessian positive definite |

**Where it can fail** (misspecification arms, one replicate each):

| simulated world | `x_s` rank | visibility rank | **`S0_c` rank** | noise |
|---|---|---|---|---|
| control, truth representable | 0.82 | 1.00 | **0.93** | 1.27x |
| each country's true R0 differs (sd log 0.10) | 0.75 | 1.00 | **0.01** | 1.18x |
| a second wave the SIR cannot make | 0.93 | 0.96 | **0.97** | 1.27x |

The country ranking collapses when transmissibility genuinely differs between countries, exactly as
the composite reading predicts; the season effects survive both violations; and the noise budget sees
neither. Same picture as under the previous model, now as the design's stated caveat rather than a
hidden one.

**The learned season pattern.** `x_s` spans -0.33 (2015/2016) to +0.55 (2025/2026) on the logit scale
-- a typical country's `S0` from 0.82 to 0.92 -- with sd 0.25 against sd 0.43 for the visibility
effect. The two are **negatively** related across seasons (correlation of the point estimates -0.57):
the seasons that spread most easily are the ones in which the smallest share of infections became a
positive consultation. Since the same-season posterior correlation is only 0.27, this is a pattern in
the estimates and not a degeneracy -- but 2025/2026, which carries the largest `x_s`, rests on 7 of 12
countries on truncated grids, so read that endpoint with its sample size.

## Fixed, not fitted

`R0 = 1.5` everywhere (see the section above); `gamma = 1/3.6` per day; `ve_inf = 0.25`, `ve_ili_cond_inf = 0.20`, `ve_spread = 0.20`; vaccination
pulse on season day 62 into 65+ only, at the reported national coverage; the contact matrix, rescaled
to spectral radius one, and re-rescaled after the elderly susceptibility re-weights it, so `R0_s` is
independent of `sigma_eld` (verified); no waning, no ageing, no importation
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

**Measured under the PREVIOUS model (`R0_s` fitted per season, `S0_c` per country) on the
86-country-season design**, before the positivity-encoding exclusion and before `R0` was fixed. These
are structural results about the estimator -- rank recovery, interval coverage, and
whether a known driver effect survives the two-step procedure -- and one country-season out of 86 does
not overturn them; but the numbers below have not been re-measured on the 85-cell design, and figure
16 is not regenerated until `run_joint_recovery.R` is run again.

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
silent -- a flat-lined country is obvious in figure 06 and in the per-country correlations -- so
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
