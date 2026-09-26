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

We fit weekly influenza-positive ILI consultations from **twelve EU/EEA countries over eight seasons —
85 country-seasons, 8 685 observed age-week cells — with a single age- and vaccination-structured SIR,
estimated jointly in 182 parameters**. **R0 is fixed at 1.5 from the literature, not fitted.** Three age groups (0-14, 15-64, 65+) mix by each country's
contact matrix, rescaled so its dominant eigenvalue is one, which makes a season's basic reproduction
number the realised R0 given that mixing rather than an artefact of the contact data's absolute scale.
The elderly susceptibility re-weights that matrix and the result is **rescaled again**, so `sigma_eld`
redistributes *who* gets infected without changing *how transmissible* the season is: `R0` keeps its
meaning whatever `sigma_eld` does, and `sigma_eld` is therefore identified by the age composition of
cases rather than by the size of the wave.
Each country-season is an independent epidemic, seeded on 1 August and integrated on a daily grid to a
fixed 53-week horizon, with a single vaccination pulse into the 65+ group on 1 October. **A country is
not one well-mixed population**: its cities are hit at slightly different times, so the national curve
is the local epidemic spread over start times drawn from a normal distribution with sd `tau` days --
one `tau` for every country and season (2026-09-26; a per-country `tau` was tested, see below).

**What varies where is the model's central design**, and it is forced by the data rather than chosen
for convenience. Observed ILI+ levels differ more than a hundredfold between countries and each country
keeps its rank across seasons, so the level has to be a **country** property — a reporting proportion —
not an epidemiological one. Seasons rise and fall **together** across Europe once that country level is
divided out, which is what licenses one Europe-wide season effect on each quantity. Transmissibility
and susceptibility enter a wave's rise rate as a **product** and are indistinguishable from one wave,
so one of them is pinned: `R0` is fixed and **`S0` is the single, deliberately blunt sensor of how
easily a wave spread** -- it reads susceptibility and infectivity together and does not try to split
them. That choice is an informed judgement, not a data result: a model comparison (below) showed the
data cannot tell a season or country effect on `S0` from one on `R0`. `S0` is decomposed additively on the logit scale into a **country level** `S0_c`
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
contracts against its prior by 0.69-0.95 (the spatial spread, which at two days barely acts, by 0.47),
and the season effect on susceptibility by 0.93, so nothing
reported is a restated assumption. The question the design must answer from shape alone -- *a bigger
season: more susceptible, or more visible?* -- it does answer: the posterior correlation between the
two season effects is only 0.27 (max 0.31 over seasons). Two limits are measured rather than asserted:
the deterministic mean needs about 2.6x more observation noise than the data's own week-to-week
scatter can explain, so some of what the model calls measurement error is really misfit; and the
country ranking of `S0_c` is **not a ranking of immunity**: on simulated data where countries' true
R0 differed by 10%, the ranking of the pure susceptibility component fell from 0.93 to 0.01, because
with R0 pinned a real transmissibility difference lands in `S0_c` -- which is exactly what a blunt
sensor is meant to do, and why it must be read as "how easily influenza spread here". A **model comparison**
asks where the variation prefers to live -- in the susceptible pool or in transmissibility -- by
moving the season effect, the country effect, or both from `S0` onto `R0` in the same 181-parameter
layout. The raw log-likelihood favours `R0` by 10-12 nats, but that number counts autocorrelated
weeks as if independent. With the **wave as the unit**, the preference is not there: 37 of 85 waves
favour `R0` (sign test p = 0.28), every cluster-bootstrap interval on the total spans zero, and the
raw gap is carried by Croatia and 2025/2026 alone -- the cells where this model's `S0` presses its
ceiling of 1 at `R0 = 1.5`. Pinning `R0 = 1.7` inside this model recovers 95% of that gap with the
country ranking preserved (rank correlation 0.99). **The data do not distinguish the two mechanisms**, so the choice is imposed
by judgement (owner, 2026-09-26): R0 fixed, `S0` the blunt sensor. The comparison's by-product is that
at 1.5 the sensor saturates (`S0` near 1) in 8 of 85 cells. And the **spatial spread** asks whether national waves are
fatter than one well-mixed epidemic makes them: a common spread of a week or more is rejected (fitted
about two days, no gain; the estimate runs low, so a few days cannot be excluded), and a per-country
spread, though it fits better in
total, is carried by a few waves, does not follow country size, and trades against `S0` within each
country -- so the working model keeps one shared tau, and a per-country tau would need outside
information on the timing of regional peaks to be identifiable. Across all three quantities that make
a wave fat, **R0 and `S0` are exactly one number; `S0` and tau are separable only through the low weeks
at either end of a wave, and in practice only for spreads of two weeks or more.**

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
| `tau` | nothing: one global value (per country in the tested variant, 12) | 1 | spatial spread, days: the sd of the start times of a country's local epidemics |

So the two quantities the project reads are each a **two-way additive decomposition**: on the logit
scale `logit S0_{c,s} = S0_c + x_s`, and on the log scale `log c_{c,s} = c_c + delta_s`. The country
levels carry the mean (their priors are centred on the plausible level, not on zero); the season
effects are centred on zero and constrained to average zero, which is what makes the level and the
effects separately identified rather than sliding against each other.

Total **182**, of which 16 are shared across countries and 166 are local to one country (with `tau`
by country: 193, of which 15 shared). The likelihood is therefore separable given the shared block,
which is what the fitting strategy exploits.
(Without the provisional positivity-encoding exclusion it is 86 country-seasons and 184 parameters —
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
influenza spread here / this season", not as a pure immunity measure.** This is the design, adopted by
judgement after the data could not separate the two (next-but-one section): `S0` is a blunt sensor of
susceptibility and infectivity together, and no parameter in the model carries R0.

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
| the data decide, not the priors | contraction 0.95 (`S0_c`), 0.93 (`x_s`), 0.94 (`c_c`), 0.90 (`delta_s`), 0.69 (`sigma_eld`); the spatial spread 0.47, because at two days it barely acts (its profile is the better description, next section) |
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

## Spatial spread, and what makes a wave fat: R0, S0 and tau (2026-09-26)

**The extension.** A country is not one well-mixed population: its cities are hit at slightly
different times. The model keeps one local epidemic per country-season and treats the country's many
local epidemics as copies of it whose start times are spread normally, sd `tau` days, around the
modelled one. The national incidence is the local incidence convolved with that normal kernel -- on
the daily grid, before the weekly aggregation, each day weighted by the normal mass in its bin and the
kernel cut at 7 sd. One parameter, `log tau ~ N(log 7 days, 0.7)`, shared by every country and
season; `tau_by_country = TRUE` gives each country its own instead. Figure 03 has a panel for it.

Three properties are exact, and the test suite holds the kernel to them: the **total** is unchanged
(so is the attack rate), the **mean timing** is unchanged, and every **exponential rate** is unchanged
-- the rise and the decline -- because a convolved exponential is the same exponential times a
constant. Only the peak changes: the curve's variance grows by exactly `tau^2 + 1/12` day^2. So tau
cannot make a wave faster, bigger or earlier; it can only round and widen it. One side effect: in the
exponential phase the spread curve runs `r tau^2 / 2` days ahead of the local one (the earliest cities
lead), which the seed absorbs.

**Checked, not assumed.** Below 0.05 days the untouched original code path runs. Against the
pre-extension engine compiled from git, the fitted means are **bit-identical in all 85
country-seasons** and the working fit's log-likelihood is identical to 0.000 nats; a refit from
scratch with tau pinned at zero lands on the old optimum (-50798.60 against -50798.59). The C++ and the
base-R reference compute the spread two different ways -- cumulative incidence at week boundaries with
erfc weights, and direct daily convolution with pnorm weights -- and agree to 1e-14 relative from
tau = 1e-13 to 119 days.

**One tau for every country: a few days at most.** Fitted tau 2.1 days from two starts; against no
spread it gains -0.08 nats (wave-bootstrap interval [-1.5, +1.2]). The profile, tau held and the other
181 parameters refitted (`run_tau_analysis.R`):

| tau held (days) | 0 | 1 | 2.1 | 4 | 7 | 10 | 13.9 | 20.5 |
|---|---|---|---|---|---|---|---|---|
| log-likelihood vs the best | -0.02 | 0 | -0.10 | -0.73 | -6.4 | -30 | -107 | -332 |

A common spread of a week or more is rejected. The estimate itself runs low (recovery, below), so the
fitted 2.1 days means "a few days at most", not "none". Not the S0 ceiling: with R0 pinned at 1.7,
where no `S0` exceeds 0.90, tau again settles at 2.3 days (+0.32 nats over no spread).

**One tau per country: real, but not geography.** Five starts reach the same optimum to 0.00 nats (a
sixth, warm-started from the shared fit, stalled 11 nats worse, which is why the runner uses five).
It fits **+52.7 nats better** than one shared tau with 11 more parameters (+35 after the
autocorrelation deflation). With the wave as the unit the evidence is weaker than the total suggests:
49 of 85 waves favour it (sign test p = 0.19, Wilcoxon p = 0.17), and the cluster bootstrap on the
total excludes zero only when whole countries are resampled ([+5, +109]; waves [-3, +124], seasons
[-4, +113]). Each country's own profile -- its tau held, its block refitted, the shared block fixed:

| country | nats lost if its tau is forced to 0 | best tau (days) | its `S0` at tau 0 -> at the best tau |
|---|---|---|---|
| PL | 20.2 | ~30 | 0.82 -> 0.88 |
| EE | 14.7 | ~30, barely bounded above | 0.80 -> 0.85 |
| NO | 8.8 | 21 | 0.85 -> 0.88 |
| FR | 6.0 | 14 | 0.90 -> 0.92 |
| DK | 5.8 | 14 | 0.88 -> 0.89 |
| IT | 3.1 | 14 | 0.88 -> 0.90 |
| CZ, BE, ES, HR, IE, NL | ~0 | 0-2 | unchanged |

Six countries want spread and six do not, and **country size does not predict which** (rank
correlation of the fitted tau with population 0.03): Spain and the Netherlands want none, Estonia the
most. Estonia's weekly series are the noisiest in the panel (dispersion 0.25) and a 36-day spread lays
a low, broad mean through the scatter; Poland is the plausible case -- without spread the model
overshoots its peaks up to twofold. So a per-country tau is a **blunt sensor of how fat a country's
waves are** -- spatial asynchrony, but equally subtype mixing, aggregation over a sentinel network
and plain noise -- not a map of its cities.

**What makes a wave fat, and what the data can tell apart.** Three quantities set a wave's shape: R0
and `S0` through the local epidemic, tau through its spread; reporting sets its height and the seed
its timing.

1. **R0 and `S0` are exactly one number.** Every term of the age x vaccination system and of its daily
   Euler step is homogeneous in the compartments, so incidence(t; R0, S0, I0) = S0 x incidence(t;
   R0 S0 at full susceptibility, I0/S0). Hence (R0, S0, I0, c) -> (k R0, S0/k, I0/k, k c) leaves every
   expected count unchanged -- verified to 7e-14 in all 85 waves (1e-9 with spread). No amount of data
   separates them, which is why fixing R0 had to be a judgement. The only asymmetry is `S0 <= 1`: the
   pin's value reaches the data through the ceiling alone (8 of 85 cells above 0.95 at 1.5), plus, weakly, the
   logit-additive form of the season and country effects, which a rescaling of `S0` does not preserve.
2. **`S0` and tau are separable in principle, through the tails.** `S0` moves the rise rate, the decline
   rate and the width together; tau moves the width only (the exact invariances above). Figure 03
   (middle of the lower row) shows it: at the same peak and width, a less susceptible wave rises AND
   decays more slowly, and by week 45 its tail is eight times higher. Two practical limits. The tails
   are the low-count weeks at either end, where the off-season baseline and the negative-binomial noise
   sit. And a weekly series' measurable rise lies partly in the approach to the peak, which spread does
   slow: in figure 18b a spread of up to a week stays within half a week of the no-spread line -- a
   small `S0` change mimics it -- and only from two weeks on does the wave leave the line. Tau is also
   one-sided: it can widen a wave, never sharpen it.
3. **What the data do with one shared tau.** The Hessian finds no correlation above 0.18 with any
   other parameter (Croatia's `S0` and seeds) and a contraction of 0.47 -- but at two days tau barely
   acts, so the profile above is the better description. Recovery on simulated data finds the ridge.
   A true shared spread of 2.1 days comes back as 2.4; of 7 days as 2.9, 3.0 and 3.0 in three
   simulations; of 14 days as 9.3, 9.5 and 10.5 -- each time with `S0` pushed down (by 0.08-0.11 logit
   at a week, 0.16-0.20 at two), while the country ranking of `S0` survives (0.94-0.99) because the
   shift is common to every country. **It is the likelihood's own ridge, not the priors**: on one of
   the 7-day data sets, holding tau at the true 7 days fits 0.64 nats WORSE than holding it at 3, and a
   threefold weaker `S0` prior leaves the estimate at 2.95 days. Up to about a week, more spread with
   higher `S0` or less of both fit almost equally well, and the maximum sits low, as a variance
   estimated from noisy data does. Read back on the real data: the profile falls 6.4 nats from the
   best tau to a week, ten times more than on data where a week is the truth, so a common week of
   spread is rejected -- but the fitted 2.1 days is compatible with a true common spread of a few days.
4. **What the data do with one tau per country.** The ridge becomes per country. Where a country uses
   spread plausibly, its tau and its `S0` are strongly correlated in the posterior -- 0.85 in France,
   0.90 in Norway, 0.70 in Italy, 0.72 in Poland, with that country's seeds at -0.7 to -0.9: more
   spread, higher `S0`, earlier seeds. Estonia's tau is pinned tightly (posterior sd 0.03 on the log
   scale) but by its noise, not by the ridge (correlation with its `S0` 0.12); the six countries near
   zero are held by the prior (contraction 0.05-0.40). In recovery from a world where big countries
   spread 14 days and small ones 3, the estimates separate in the right direction (medians 9.6 against
   4.6 days, rank correlation 0.77) but shrink towards each other, and the **country ranking of `S0`
   comes back at 0.69 instead of 0.94-0.99**, the season effects at 0.68 instead of 0.82-0.93. The
   wave-level test CAN see a real country difference: on those simulated data 57 of 85 waves favour
   the per-country model (sign test p = 0.002, Wilcoxon p = 0.0001) and every bootstrap interval
   excludes zero, for a modest +12 nats. On the real data the total is four times larger but only 49
   waves agree -- a large gain carried by a few waves, the signature of absorbing wave-specific misfit
   rather than of a country-level spread.
5. **What moves when tau is by country, and what does not.** On the real data the country `S0` ranking
   moves (rank correlation 0.80 against no spread; Estonia from lowest to fourth, the Netherlands from
   fourth to lowest), so does the 2025/2026 season effect (0.55 -> 1.01 logit) and the number of `S0`
   at the ceiling (8 -> 10). The pattern of the season effects on `S0` does not (0.98), nor do season
   visibility and reporting levels (0.99). And the pin reaches tau nowhere: with R0 at 1.7 every
   country's tau moves by at most a day -- Croatia, at the `S0` ceiling under 1.5, stays at 2.5 days
   although its local wave may now be sharper -- while every `S0` scales by 1.5/1.7 = 0.88 exactly as
   point 1 requires (Denmark 0.887 -> 0.782), Croatia excepted because the ceiling released it
   (0.972 -> 0.867).

**Decision.** The working model keeps **one shared tau** (`tau_by_country = FALSE`, the default):
the data put it at about two days, so every conclusion of the model without spread stands, and the
parameter now reports a result -- no common spread of a week or more. **Tau by country is not
adopted.** It fits better in total but not wave by wave, where a real country difference would show;
its values do not follow country size and in Estonia fit noise; and within each country it trades
against `S0` -- in the simulation where it is the true model, the country ranking of `S0` came back at
0.69, against 0.94-0.99 whenever tau is shared -- which costs the quantity the project reads. It stays
one setting away as a tested hypothesis (`run_tau_analysis.R`, figure 18). **If spatial spread per
country is to be modelled, it needs outside information to pin it**: the observed dispersion of
regional peak times within each country, from sub-national surveillance or published estimates, as a
prior on each tau. That would take it off the ridge instead of letting it float there, and it is what
more of the same national data cannot do.

## Where does the variation live: in `S0` or in `R0`? (model comparison, 2026-09-25)

**Outcome (owner, 2026-09-26): the data cannot decide, so judgement does.** R0 stays fixed and `S0`
is a blunt sensor of susceptibility and infectivity. The sensing switch that made this comparison
possible has been removed from the code to keep the model simple; the comparison is reproducible from
commit `7b04681` (`run_model_comparison.R`, `compare_inference.R`, `plot_model_comparison.R`). What
follows is the evidence behind the decision.

**The question.** The working model senses season and country with the susceptible pool:
`logit S0_{c,s} = S0_c + x_s`, `R0 = 1.5`. The alternative senses them with transmissibility:
`log R0_{c,s} = R0_c + r_s`, `S0 = 0.75`. The two are the same layout with one anchor and one fitted
family swapped -- 181 parameters either way -- so each effect can be placed on either quantity
independently, giving four models: `S0/S0` (working), `R0/S0`, `S0/R0`, `R0/R0` (season / country).
Same data, same likelihood, not nested: a likelihood comparison, not a test. What the likelihood
scores is shape: an effect on `S0` changes how fast a wave rises *and* how many are left to infect
(a deeper pool decays differently); an effect on `R0` changes the speed only. Figure 17
(`output/joint_model/17_model_comparison.png`, produced at commit `7b04681`).

**The raw numbers, and why they overstate.**

| model | season on | country on | log-likelihood | vs working | AIC vs working |
|---|---|---|---|---|---|
| `S0/S0` | S0 | S0 | -50798.59 | 0 | 0 |
| `R0/S0` | R0 | S0 | -50788.43 | +10.16 | -20.3 |
| `S0/R0` | S0 | R0 | -50786.95 | +11.64 | -23.3 |
| `R0/R0` | R0 | R0 | -50786.92 | +11.67 | -23.3 |

Either lever on `R0` gains 10-12 nats; both together gain nothing more (interaction -10.1 nats), so
this is one signal, not two. But nats and AIC count each week as an independent observation, and the
weeks inside a wave are not: the residual lag-1 autocorrelation within waves has median 0.21
(quartiles -0.01 and 0.42), an AR(1) effective-sample-size factor of 0.65, so the honest size of the
gap is nearer 7 nats -- and the unit of inference has to be the wave.

**With the wave as the unit** (85 country-seasons; per-wave paired log-likelihood differences; sign
test, Wilcoxon signed-rank, and a cluster bootstrap of the total resampling waves, countries or seasons):

| contrast (`R0` variant minus working) | waves favouring `R0` | sign test | Wilcoxon | total | 95% CI, waves | countries | seasons |
|---|---|---|---|---|---|---|---|
| season effect on `R0` (`R0/S0`) | 37 of 85 | p = 0.28 | p = 0.61 | +10.2 | [-17.8, +44.5] | [-26.4, +59.6] | [-17.9, +41.9] |
| country effect on `R0` (`S0/R0`) | 39 of 85 | p = 0.52 | p = 0.67 | +11.6 | [-16.6, +46.3] | [-26.5, +64.4] | [-16.8, +42.9] |
| both on `R0` (`R0/R0`) | 38 of 85 | p = 0.39 | p = 0.65 | +11.7 | [-16.6, +46.4] | [-26.7, +64.5] | [-16.8, +43.0] |

The *typical* wave prefers `S0` (median per-wave difference -0.04 to -0.06 nats); the total is
positive because a few waves prefer `R0` by a lot. Croatia carries +18.9 of the +10.2 and 2025/2026
+12.4 (Croatia 2025/2026 alone +12.0); every other cell together nets -9.1, against `R0`.

**Why those cells: the ceiling.** With `R0 = 1.5`, the fastest rise the working model can produce is
`gamma (1.5 - 1)` -- at `S0 = 1`. It presses that ceiling: 8 of 85 fitted `S0_{c,s}` exceed 0.95
(maximum 0.985; Croatia's country level 0.974). The `R0` variants have no ceiling. Refitting the
working model with a higher pin and nothing else changed:

| pinned `R0` | vs `R0 = 1.5` | max `S0_{c,s}` | cells above 0.95 | Croatia level | gain in Croatia | in 2025/2026 | all other cells |
|---|---|---|---|---|---|---|---|
| 1.5 | 0 | 0.985 | 8 | 0.974 | | | |
| 1.6 | +7.00 | 0.944 | 0 | 0.917 | +8.4 | +6.4 | -2.1 |
| 1.7 | +9.64 | 0.900 | 0 | 0.864 | +12.3 | +9.1 | -3.6 |

`R0 = 1.7` recovers 95% of what the season-on-`R0` model gains, in the same cells, with no effect
on `R0` at all, and the country ranking of `S0_c` is preserved (Spearman 0.993 between the 1.5 and
1.7 fits). What the `R0` variants "know" is only that Croatia and 2025/2026 rose faster than a pool
of at most 100% susceptibles can at `R0 = 1.5`.

**Verdict.** The data do not distinguish sensing with `S0` from sensing with `R0`: a rise-rate
difference can be booked to either, and the shape signal that separates them is worth a few nats
spread over 85 waves, pointing slightly to `S0` in the typical wave and strongly to `R0` only where
`S0` has run out of room. The model stays `S0/S0` by judgement. What the comparison did
establish is that **the sensor saturates at `R0 = 1.5`**: `S0` presses 1 in 8 of 85 cells, and a pin
of 1.7 frees every cell (no `S0_{c,s}` above 0.90) at the cost of every `S0_c` moving down (Croatia
0.97 to 0.86), ranking unchanged. The pin stays 1.5 until the owner decides otherwise; the blunt
reading of `S0` is the same at either value.

**What the comparison also found.** The first `R0 = 1.6` refit reported a log-likelihood of +3 000 000:
Spain's dispersion had run to `log phi = 44.7`, where `lgamma(y + phi) - lgamma(phi)` is catastrophic
cancellation and comes out hugely positive. Both the C++ and its base-R mirror use that formula, so
the identity test was blind to it. Closed by rejecting `phi > 1e8` in both (the negative binomial is
Poisson to eight digits there, and every real fit has `phi` between 0.15 and 1.8), with a test
against R's own `dnbinom` as an independent third implementation. No previous result was in that
region; the refit above is the valid one.

## Fixed, not fitted

`R0 = 1.5` everywhere (see the two sections above); `gamma = 1/3.6` per day; `ve_inf = 0.25`, `ve_ili_cond_inf = 0.20`, `ve_spread = 0.20`; vaccination
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
