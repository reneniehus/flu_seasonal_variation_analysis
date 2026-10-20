# The joint EU/EEA influenza model

**How to read this file.** Sections 1-7 describe the model as it stands -- R0 pinned at 2.0, spatial
spread fixed at none, children's susceptibility raised by a modifier anchored on cohort evidence, no
susceptibility factor for the elderly, 181 parameters, fitted 2026-09-27 -- and what it has learned, with every number taken from that fit
(`output/joint_model/joint_fit.rds`). Section 8 is a guide to the figures. The
appendices are the evidence trail: the investigations that shaped the design, each dated, with its
numbers as measured at the time. `documentation/decisions.md` records every decision in order, and
`documentation/to_confirm_with_surveillance.md` the six things we have inferred about the source data
rather than confirmed. For reading on paper, `documentation/joint_model_report.pdf` tells the same
story in 16 pages (built by `code/07_joint_model/report/`).

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

## 1. In one page

We fit weekly influenza-positive ILI consultations (ILI+) from **twelve EU/EEA countries over eight
seasons -- 85 country-seasons, 8 685 observed age-week cells -- with one age- and
vaccination-structured SIR per country-season, estimated jointly in 181 parameters.** The model is
best read as a **structured sensor of wave shape**. It turns each country's weekly curves into two
things: how easily influenza spread, a susceptibility-like index `S0` that varies by country and by
season, and how visible that spread was, a reporting level that varies by country and by season. The
SIR ties each wave's rise, peak and decline to one number, which is what lets the two be told apart.

**Three features of the data force the design** (figure 02):

1. **The level is a country property.** Peak ILI+ differs about 150-fold between countries and each
   country keeps its rank across seasons: that is how much a surveillance system sees, not how many
   people were infected. So each country has its own reporting level.
2. **Seasons move together across Europe.** Once each country's level is divided out, a big season is
   big almost everywhere. So each season has ONE effect on `S0` and ONE on visibility, shared by all
   countries -- which is also what makes the joint fit worth doing.
3. **Waves arrive at different times.** Peaks span 14 weeks across the panel. So every country-season
   has its own seed, the model's only handle on timing.

**Two quantities are pinned by judgement, not fitted**, because the data cannot separate them from `S0`:

- **R0 = 2.0.** R0 and `S0` are exactly one number to the data: multiplying R0 by any factor while
  dividing `S0` and the seed by it and multiplying reporting by it reproduces every expected count. The data fix the
  product -- a reproduction number at season start of 1.28 in the median country-season, where the
  literature puts that of seasonal influenza -- and the pin chooses only how it splits. Its value sets the absolute scale of
  infections (attack rates scale as 1/R0) and no relative result; 2.0 keeps every `S0` well below its
  ceiling of 1.
- **Spatial spread = none.** A country's cities are hit at slightly different times, which fattens
  the national wave. A fitted spread traded against `S0` and cost its country ranking its
  identifiability, so it is fixed and `S0` absorbs it.

**One more rests on outside evidence: children start each season more susceptible.** They carry the
least immunity from past seasons, and cohorts that test everyone find them infected about 1.7 times as
often as adults (PHIRST). One global modifier raises children's `S0` on the logit scale and can only
raise it. The surveillance data cannot tell children's susceptibility from how readily children are
seen, so the modifier's prior is anchored on the cohort evidence and the data pull it partway back:
children's `S0` is 0.79 against 0.59 for adults, and they are infected 1.42 times as often. The elderly
get no factor of their own: the same cohorts find them infected no more often than adults -- ageing
raises severity, not infection -- so per contact they are as susceptible as adults.

So **`S0` is a blunt sensor**: susceptibility, infectivity (a more transmissible strain, denser mixing,
an older population) and the overlay of a country's local waves all land in it. Read `S0_c` and the
season effects as "how easily influenza spread here" and "this season", never as immunity alone.

**What it learned** (figures 09-14):

- **Seasons differ little in how easily they spread and a lot in how visible they were.** A typical
  country's adult R0 x `S0` was 1.12 in 2015/2016, the lowest, and 1.30 in 2025/2026, the highest; the
  other six seasons lie between 1.17 and 1.20, closer than the data can order. Visibility per infection ran
  from 0.57 times normal (2025/2026; 0.58 in 2023/2024) to 1.84 (2017/2018), about ten times the spread of the
  susceptibility effect on a common log scale. The two move against each other (correlation -0.61):
  the seasons that spread most easily turned the fewest infections into consultations.
- **Countries.** Croatia's waves spread most easily (adult R0 x `S0` 1.40), Ireland's and Estonia's
  least (1.11, 1.13). The ranking is sharp -- intervals of about +/-0.01 to 0.03 on `S0` -- but it
  assumes equal transmissibility and equal spatial structure across countries, which is exactly what
  the pins impose.
- **Reporting.** The share of adult infections that becomes an ILI+ consultation differs about
  twentyfold between countries, Italy and Poland highest, Croatia and Ireland lowest. A child is seen
  about 1.24 times as readily as an adult per infection -- most of children's excess cases are more
  infection, not more reporting -- and an elderly person about 1.8 times: the elderly's excess cases
  are visibility, which fits their more severe illness.
- **Who gets infected.** Children most (34% in the median country-season), adults 22%, the elderly
  6.5%: children 1.42 times adults, the elderly 0.29 times, against 1.7 and 0.8 in PHIRST. With no
  elderly factor the contact matrices alone set the elderly's share, well below the cohort's -- whose
  elderly live in multigenerational South African households -- and the data would raise it slightly
  (section 6). The absolute level, 20% overall, is the pin's choice; the pattern is the data's and the
  cohorts'.

**How far to trust it** (sections 5 and 6). Every parameter family is determined by the data rather
than its prior (contraction 0.78-0.97), except the children's modifier, which is anchored on outside
evidence by design (0.23). The two questions the design must answer from shape alone --
*a bigger season: more susceptible or more visible?* and *a big country: more infections or more
seen?* -- are answered: the posterior correlations are 0.23 and 0.10. On data simulated from the model onto the real design and refitted
from scratch (re-run 2026-09-27 on this model), the country ranking of `S0` comes back at 0.97, the
season effects at 0.97 and visibility at 0.99, the intervals cover 93% of the time, and a known driver
effect on the season parameters is recovered without bias. The children's modifier comes back where its
prior pulls it, not where the truth is, and the level of `S0` shifts with it (section 6). Three
limits are measured, not assumed: the deterministic curve needs 2.6 times more noise than the data's
own week-to-week scatter, so part of what it calls measurement error is misfit; the country ranking of
`S0` does not survive real transmissibility differences between countries (it falls to 0.27 on
simulated data where they differ by 10%); how much more susceptible children are rests on the
cohort evidence, not on these data; and without an elderly factor the contact matrices set the
elderly's share of infection, which the data would raise slightly (10.9 nats; 58 of 85 waves). One pattern is unexplained: the country ranking of `S0` follows
data quality (rank correlation 0.52 with the dispersion; 0.75 before the children's modifier), and
recovery shows the estimator does not produce it.

## 2. The data the model sees

- **Source.** Weekly ILI+ per 100 000 by age group (0-14, 15-64, 65+): RespiCompass for the seasons up
  to 2023/2024, ERVISS (ILI rate x sentinel positivity; non-sentinel positivity for Croatia) from
  2024/2025, stitched per week into one panel (`stitch_iliplus.R`, figure 01). Rates are converted to
  counts with each age group's own population, a scale device the reporting level absorbs; the
  negative binomial then works on counts.
- **Design.** 12 countries (DK, EE, ES, FR, NO, BE, CZ, IE, IT, PL, HR, NL), 8 seasons (2014/2015 to
  2018/2019 and 2023/2024 to 2025/2026, the COVID seasons omitted), 85 country-seasons on grids of 33
  to 53 weeks from 1 August, of which 16 to 53 weeks carry data -- 8 685 observed age-week cells, 31.5%
  of them exactly zero. A country
  needs five age-complete seasons to enter (Spain has five, Norway and Czechia six).
- **Excluded, provisionally.** Czechia 2024/2025, because of how ERVISS encodes weeks without a
  detections row (box above); everything else is kept.
- **Uneven support.** 2025/2026 rests on 7 of the 12 countries (DK, EE, FR, BE, IE, PL, HR), on grids of
  34-41 weeks, with 774 observed cells against 1 497 for 2023/2024 -- and it is the season with the
  largest effect. Read it with its sample size (`jm_summary_season` reports it beside every number).
- **Mixing.** Each country has its own Prem et al. contact matrix except Norway, which uses the EU
  average; the age reporting offsets are the parameters most sensitive to that.
- **A deliberate selection.** At the design's own inclusion rule Iceland, Malta and Austria would also
  qualify and are not offered to the model. Since the season effects are shared across the countries
  in the design, which countries are in it is substantive.
- **Audited** (appendix C): every observation matrix rebuilds exactly from the raw streams by
  independent code; age bands and populations agree to the person.

## 3. The model

**Dynamics.** Each country-season is an independent epidemic: three age groups x vaccinated or not,
SIR, integrated on a daily grid (forward Euler) from a seed on 1 August to a fixed 53-week horizon,
with one vaccination pulse into the 65+ group on 1 October at the reported national coverage. The
country's contact matrix is rescaled to spectral radius one, so the realised reproduction number at
full susceptibility is exactly the pinned R0 in every country. Infection per contact is the same at
every age; only children start a season with a higher `S0`.

**Observation.** Weekly counts are negative-binomial around the model mean, which is (infections in
that age group and week) x (reporting level of the country x visibility of the season x the age
offset) + an off-season floor per data source, with a dispersion per country.

**What varies where.**

| parameter | varies | count | meaning |
|---|---|---|---|
| `S0_c` | between countries | 12 | the country's `S0` at the average season, logit scale |
| `x_s` | between seasons, shared by all countries, averaging zero | 7 free | season effect on logit `S0` |
| `c_c` | between countries | 12 | the country's reporting level at the average season, log scale |
| `delta_s` | between seasons, shared, averaging one | 7 free | season effect on visibility, log scale |
| `kappa` | nothing: one number, only positive | 1 | children's `S0` modifier: `logit S0_young = logit S0 + kappa` |
| `off_c,young`, `off_c,eld` | between countries | 24 | how readily children and the elderly are seen, relative to adults |
| `phi_c` | between countries | 12 | negative-binomial dispersion |
| `b_c,src` | between countries and data sources | 21 | off-season floor per source present |
| `I0_c,s` | freely between country-seasons | 85 | the seed, i.e. when the wave arrives |

So `logit S0_{c,s} = S0_c + x_s` and `log c_{c,s} = c_c + delta_s`: two two-way additive designs on 85
cells. The country levels carry the mean; the season effects average zero, which is what separates
them. Children's `S0` sits `kappa` above the adults' on the logit scale in every country and season.
181 parameters, 15 shared across countries and 166 local to one, so the likelihood is separable
given the shared block and is fitted by block coordinate descent (appendix I).

**Fixed, not fitted.** R0 = 2.0 everywhere and spatial spread 0 (section 1; appendices E and F);
infectious period 3.6 days; vaccine effects on infection 0.25, on ILI given infection 0.20, on onward
spread 0.20; contact matrices as above; no waning, no ageing, no importation after the seed, no latent
period, one strain per season. Priors are weak and centred on plausible values; the `S0` prior is
placed on R0 x `S0` (centre 1.125) so the pin cannot reach the fit through it. The one informative
prior is on `kappa`: log `kappa` ~ N(log 1.5, 0.2), centred where the refitted model reproduces
PHIRST's children-to-adult infection ratio of about 1.7 (appendix L).

**What each fitted quantity means now.**

- `S0_{c,s}`: how easily influenza spread in that country and season, as a susceptible fraction at
  R0 = 2. It carries susceptibility, strain transmissibility, a country's demography and mixing (the
  contact matrix is rescaled, so real differences in overall transmissibility land here), and the
  overlay of the country's local waves. Its ranking is the reading; its level is the pin's.
- `c_c` and `delta_s`: how many ILI+ consultations one infection produces in that country, and how much
  more or less in that season (how symptomatic the strain was, how many consulted, how many were
  tested). Their product is all the data see; the average-one constraint on `delta_s` separates them.
- `kappa`: how much less immunity children carry into a season than adults, as a shift on the logit
  scale, one number for Europe. Biology, anchored on cohort evidence; the surveillance data alone
  cannot separate it from children's reporting.
- The elderly have no factor of their own: cohorts find them infected no more often than adults --
  ageing raises the risk of severe outcomes, not of infection (appendix L) -- so per contact they are as
  susceptible as adults, and their infections follow from their contacts and their vaccination.
- The age offsets: surveillance -- how readily each age group is seen, per country.

## 4. What it learned (the working fit, R0 = 2.0, 2026-09-27)

**Seasons** (figures 09 and 10; intervals are 95% from the curvature of the fitted surface):

| season | countries | adult `S0`, typical country | adult R0 x `S0`, typical | visibility x normal (95%) | median attack rate |
|---|---|---|---|---|---|
| 2014/2015 | 11 | 0.588 | 1.18 | 1.40 (1.28-1.54) | 21% |
| 2015/2016 | 12 | 0.562 | 1.12 | 1.16 (1.05-1.27) | 16% |
| 2016/2017 | 12 | 0.597 | 1.19 | 0.93 (0.85-1.02) | 21% |
| 2017/2018 | 12 | 0.585 | 1.17 | 1.84 (1.68-2.02) | 20% |
| 2018/2019 | 11 | 0.593 | 1.19 | 1.09 (0.99-1.21) | 22% |
| 2023/2024 | 11 | 0.592 | 1.18 | 0.58 (0.52-0.65) | 19% |
| 2024/2025 | 9 | 0.598 | 1.20 | 0.99 (0.88-1.11) | 20% |
| 2025/2026 | 7 | 0.648 | 1.30 | 0.57 (0.50-0.65) | 29% |

The season effect on logit `S0` has sd 0.101; on the log scale of R0 x `S0` that is about 0.039,
against 0.40 for visibility. Only two seasons stand apart in how easily they spread: 2015/2016 below,
2025/2026 above; the other six are within the data's resolution of each other, so their order swaps
under any small change (it does between pins). Visibility is where the seasons differ: 2017/2018 and
2014/2015 turned many more infections into consultations, 2023/2024 and 2025/2026 far fewer. The
negative correlation between the two (-0.61) is in the estimates, not forced by the model (the
same-season posterior correlation is 0.23), and 2025/2026 carries both extremes on the thinnest data.
Children's `S0` is higher in every season by the same logit shift: 0.79 in the median country-season
against 0.59 for adults. With it, the reproduction number at season start is 1.28 in the median
country-season, as without the modifier: the data fix how fast the epidemic grows, and the modifier
only moves who carries it. The season and visibility effects are unchanged by it (correlation 0.999
with the fit without it); the country ranking of `S0` moves somewhat (0.90), most for Ireland, where
the modifier also costs the most fit (appendix L). Removing the elderly factor moved none of these:
season effects, visibility and the country ranking correlate 1.000 with the fit that had it.

**Countries** (figures 11 to 13):

| country | adult `S0` (95%) | R0 x `S0` | adult reporting | child x adult | 65+ x adult | dispersion | noise excess | median attack |
|---|---|---|---|---|---|---|---|---|
| HR | 0.701 (0.690-0.713) | 1.40 | 1.0% | 2.23 | 1.08 | 1.79 | 2.5x | 36% |
| ES | 0.649 (0.636-0.661) | 1.30 | 2.6% | 1.55 | 1.78 | 1.33 | 2.0x | 27% |
| FR | 0.624 (0.608-0.639) | 1.25 | 7.5% | 1.16 | 1.81 | 0.60 | 2.6x | 25% |
| CZ | 0.623 (0.611-0.635) | 1.25 | 2.1% | 2.24 | 1.46 | 0.37 | 2.5x | 24% |
| DK | 0.605 (0.587-0.623) | 1.21 | 3.1% | 1.33 | 2.78 | 0.40 | 2.4x | 22% |
| BE | 0.600 (0.582-0.617) | 1.20 | 14.2% | 0.75 | 1.85 | 0.29 | 3.9x | 21% |
| IT | 0.592 (0.567-0.616) | 1.18 | 24.5% | 1.07 | 1.90 | 0.98 | 3.5x | 19% |
| NO | 0.587 (0.573-0.600) | 1.17 | 6.1% | 0.24 | 1.28 | 0.65 | 2.3x | 19% |
| PL | 0.583 (0.570-0.596) | 1.17 | 18.3% | 2.73 | 1.93 | 0.15 | 3.0x | 18% |
| NL | 0.574 (0.551-0.597) | 1.15 | 3.4% | 1.08 | 8.63 | 0.18 | 2.5x | 17% |
| EE | 0.565 (0.551-0.578) | 1.13 | 6.8% | 2.77 | 1.25 | 0.24 | 3.4x | 15% |
| IE | 0.557 (0.529-0.586) | 1.11 | 1.4% | 0.34 | 3.40 | 0.75 | 3.0x | 19% |

**Age.** Children are infected most: 34% in the median country-season, against 22% for adults and 6.5%
for the elderly -- children 1.42 times adults (interquartile 1.30-1.60), the elderly 0.29 times;
PHIRST found about 1.7 and 0.8. `kappa` = 0.94 (95% 0.70-1.28), pulled down by the data from its prior
centre of 1.5. The reporting offsets say a child is seen 1.24 times as readily as an adult per
infection (median over countries) and an elderly person 1.83 times, from 1.08 in Croatia to 8.6 in the
Netherlands: with infection per contact the same at every age, the elderly's excess of cases is
visibility, which fits their more severe illness (appendix L).

**Fit.** Observed and fitted weekly counts correlate at a median of 0.91 per country-season (interquartile
0.85-0.95). The worst five are Estonia 2014/2015 (0.35), Estonia 2023/2024 (0.42), Czechia 2018/2019
(0.46), Estonia 2015/2016 (0.58) and Estonia 2016/2017 (0.64): Estonia's series are the noisiest in the
panel. The noise budget (figure 08): the fit needs 2.6 times the data's own week-to-week scatter
(median over countries, range 2.0-3.9).

## 5. What the data can and cannot tell apart

**The ties, and how each is broken.**

| pair | relation | broken by |
|---|---|---|
| R0 and `S0` | exactly one number (the whole wave depends on R0 x `S0`; `S0` also scales the height, which reporting absorbs) | the pin, R0 = 2.0 -- a judgement |
| reporting level and season visibility | exactly one number (a product) | the average-one constraint on `delta_s` |
| `S0` and spatial spread | nearly one number: separable only through the low weeks at either end of a wave, and in practice only for spreads of two weeks or more | fixing the spread at none -- a judgement |
| children's susceptibility and children's reporting | nearly one number: both set the level of children's counts; they differ only in how much the children's wave leads and drives the adults', which the data barely resolve (posterior correlation -0.31; appendix L) | the prior anchored on cohort evidence |
| susceptibility effect and visibility effect of a season | different shapes: a susceptible season rises faster and peaks earlier, a visible one is the same curve scaled | the data (figure 03) |
| country `S0` and country reporting | different shapes, the same way | the data |

**What the data do with the rest** (penalised Hessian of the working fit):

| claim | number |
|---|---|
| the data decide, not the priors | contraction 0.97 (`x_s`), 0.97 (`S0_c`), 0.94 (dispersion), 0.94 (`c_c`), 0.90 (`delta_s`), 0.83 (seeds), 0.82 (baselines), 0.78 (age offsets); the children's modifier 0.23, by design |
| more susceptible or more visible is separable | posterior correlation of a season's two effects: median 0.23, at most 0.32 |
| country level and reporting level are separable | posterior correlation: median magnitude 0.10, at most 0.56 |
| `S0` is read off the shape of the rise | Spearman correlation of the observed early growth rate with the fitted `S0_{c,s}`: 0.46 over 83 waves |
| no direction the data cannot see | no near-flat eigenvalue of the likelihood-only Hessian; the penalised Hessian is positive definite |

**Where it would fail** (misspecification arms, re-run 2026-09-27 on this model, one replicate each):
data simulated from truths the model CANNOT represent, refitted, and the reported rankings scored.

| simulated world | season `S0` rank | visibility rank | **country `S0` rank** | reporting rank | noise excess |
|---|---|---|---|---|---|
| control, truth representable | 0.89 | 1.00 | **0.99** | 1.00 | 1.25x |
| each country's true R0 differs (sd 0.10 on the log scale, i.e. 0.84-1.20) | 0.93 | 0.93 | **0.27** | 0.93 | 1.27x |
| a second wave the SIR cannot make (+60% late in every season) | 0.79 | 1.00 | **0.99** | 0.99 | 1.27x |

The country ranking collapses when transmissibility genuinely differs between countries -- with R0
pinned, that difference has nowhere to go but `S0_c`, which is what a blunt sensor does (0.07 in the
previous run, 0.27 in this one: one replicate each) -- while the season results largely survive both
violations. The noise budget sees neither, so it is not the diagnostic for
this; a per-country R0 fitted as an alternative model and compared wave by wave would be.

## 6. How far to trust it

**Does it recover a known truth?** (`run_joint_recovery.R`, re-run 2026-09-27 on this model, figure 16.)
Data are simulated from the model onto the real design -- the same countries, seasons, weeks, missing
cells and populations -- and the whole pipeline is refitted from scratch.

| study | what the truth is | what comes back |
|---|---|---|
| local, 8 replicates | the fitted optimum | country `S0` ranking 0.97 (median Spearman); season effects on `S0` Pearson 0.97 (Spearman 0.89: the six near-tied seasons swap); visibility Pearson 0.99 (Spearman 1.00); 95% intervals cover 93% (median over families; season effects 98%, visibility 93%, reporting 94%, but `S0_c` 77% and the children's modifier 4 of 8) |
| prior-space, 4 replicates | drawn from the priors | season effects Spearman 1.00 and 0.98; country `S0` ranking 0.997; coverage 89-100% in every family |
| driver, 6 replicates | a covariate moves both season effects by a known amount | effect on logit `S0`: 0.245 +/- 0.029 for a true 0.25; on log visibility: 0.359 +/- 0.068 for a true 0.35 -- unbiased |

What this says. On data like ours the estimator finds the season effects, visibility, reporting and
the country ranking, and their intervals mean what they say; over the wider prior space the rankings
hold and the intervals do. The learning layer's own procedure -- fit, then regress the fitted season
effects on a driver -- returns the driver's effect without attenuation, so its slopes can be read at
face value.

Two things are not recovered, and they are one thing. The children's modifier, simulated at 0.94, comes
back at 1.22 (1.10-1.36), pulled towards its prior centre of 1.5, and its intervals cover 4 of 8: an
anchored prior on weak data does exactly that. Every country's `S0` then comes back about 0.06 logit low
to keep the waves' growth -- the replicate's mean `S0` error and its modifier error correlate at -0.97 --
so `S0_c`'s intervals cover only 77%. With that common shift removed the country errors are 0.03 logit
and the ranking is untouched. So the modifier and the LEVEL of `S0` are set by the prior on the
modifier, as the design intends; the ranking and the season effects are the data's. The off-season
baseline is biased too (a nuisance parameter).

**The data-quality pattern is not the estimator's.** The country ranking of `S0` follows the
dispersion (0.75 on the fit without the children's modifier, 0.52 with it: the three noisiest series,
Poland, the Netherlands and Estonia, carry three of the four lowest `S0`). If noise flattened the fitted
waves, the recovery would show the noisy countries' `S0` biased low against the others; relative to the
common shift, Poland and Estonia come back 0.03 logit high and the Netherlands 0.02 low. So the association is in
the data, or in misfit the simulation does not contain -- an open question (section 7), and a reason
to read the bottom of the country ranking with care.

**The robustness record** (appendices):

- **Sensing with `S0` or with R0** (appendix G): the data cannot tell a season or country effect on
  `S0` from one on R0 -- with the wave as the unit no preference survives -- so the choice is a
  judgement.
- **Spatial spread** (appendix E): a common spread of a week or more is rejected; a per-country spread
  fits better in total but not wave by wave, does not follow country size, and costs `S0`'s country
  ranking. Fixed at none.
- **The pin** (appendix F, measured before the age changes of 2026-09-27): from R0 = 1.7 up the fit is flat
  within 2 nats and every relative result is unchanged (country ranking 0.99, season effects 0.996,
  visibility 0.998); only the absolute scale moves. The modifier's prior is calibrated at the pin
  (appendix L).
- **The season asymmetry is the data's** (re-run 2026-09-26 on the fit just before the children's
  modifier and the removal of the elderly factor, neither of which moves the season effects): widening the prior on the
  season effect on `S0` fourfold, or both season priors fourfold, leaves both spreads where they were
  (0.036 on the log scale of R0 x `S0` against 0.41 for visibility); tightening the visibility prior
  fourfold -- a deliberately hostile setting that costs 7 nats of fit -- still leaves visibility varying
  nine times as much (0.33 against 0.038). `output/joint_model/prior_stress_R0_2.0.rds`.
- **The children's modifier** (appendix L): the data alone prefer none (profile maximum at zero), but
  at the fitted value the cost is not evidence with the wave as the unit (at `kappa` = 1, 36 of 85 waves
  favour it, sign test p = 0.19, every bootstrap interval spans zero); it falls on two countries,
  Ireland and Croatia. Only large values (`kappa` of 2-3, children infected 2-2.4 times adults) are
  rejected, and only at the margin (Wilcoxon p = 0.03).
- **The elderly factor, removed** (appendix L): the data would keep it, modestly but consistently --
  10.9 raw nats, and 58 of 85 waves and 11 of 12 countries fit better with it (sign test p = 0.001,
  every bootstrap interval excludes zero). Nothing downstream moves: season effects, visibility and the
  country ranking correlate 1.000 with the fit that had it.
- **The noise gap is real** (appendix K): across 18 ways of measuring the data's own scatter the excess
  runs 2.1-3.7 times.
- **A defect found and closed on the way**: a region of the dispersion where the likelihood formula
  cancels catastrophically (appendix G); the likelihood is now also checked against R's own `dnbinom`.

**Caveats to carry into any use of the output.**

1. The country ranking of `S0` assumes equal transmissibility and equal spatial structure across
   countries; a 10% real difference in transmissibility is enough to scramble it (section 5).
2. Absolute attack rates and reporting fractions are set by the R0 pin; only their patterns are data.
3. 2025/2026 is the most extreme season on the thinnest data.
4. The deterministic curve does not follow the within-season wander of the data (the 2.6x noise gap).
5. How much more susceptible children are rests on PHIRST, a South African cohort, not on these data.
6. The bottom of the country ranking follows data quality: the three noisiest series carry three of the
   four lowest `S0`.
7. The elderly's share of infection is the contact matrices' (0.29 times adults, against 0.8 in
   PHIRST), and their counts are matched by visibility instead; the data would move a little of that
   back into infection.

## 7. Open questions for reflection

**About the data.**
- The six surveillance questions (`to_confirm_with_surveillance.md`), first among them whether an
  ERVISS week with tests and no detections row is a zero.
- Whether the twelve countries should grow to fifteen (Iceland, Malta, Austria qualify).
- 2025/2026's series stop after 34-41 weeks in the panel; its extremes may move once the rest of the
  season is added.

**About the model.**
- **`S0` and data quality.** The country ranking of `S0` follows each country's dispersion (0.52), and
  the estimator does not produce that (section 6). Either countries with noisier surveillance really
  have slower, flatter waves, or their misfit has a structure the simulation lacks. Worth a look before
  the country ranking is used.
- **Country x season.** The design gives each season one effect for all of Europe, so a country whose
  season departs from the European pattern cannot be fitted: in 2015/2016 Ireland ran far above the
  model and Croatia far below it (figure 06). That misfit goes into the dispersion and is part of the
  noise gap below. A country-specific season deviation, shrunk towards zero, is the natural extension --
  and the question the learning layer will face anyway (is a driver's effect the same everywhere?).
- **Seeding on 1 August.** Every wave is seeded on the same date and its seed size sets its timing. For
  Croatia, whose waves rise fastest, that takes seeds near 1e-11, three prior standard deviations below
  the rest: the device is stretched there. Seeding at the observed onset instead would decouple timing
  from growth.
- **The 2.6x noise gap.** Something pervasive across each season is missed -- school holidays, a second
  subtype, reporting behaviour. Candidates: process noise (the Kalman filter the pilot had), a
  second strain, or a holiday contact effect. The misspecification test says a late second wave alone
  would not show in the noise budget, so the gap is not that.
- **Children.** Their extra susceptibility now enters as `kappa`, anchored on PHIRST because the
  surveillance data cannot separate it from children's reporting. PHIRST is a South African cohort --
  a younger population, higher HIV prevalence, different contacts -- so a European anchor (Flu Watch in
  England, or paediatric serology) would be the better prior if one can be extracted.
- **The elderly.** With no factor, the contact matrices set the elderly's share of infection (0.29 times
  adults, against 0.8 in PHIRST, whose elderly live in multigenerational households) and a reporting
  offset of about 1.8 matches their counts. The data prefer a correction, modestly (appendix L). If the
  elderly's attack rate ever becomes a target, the correction belongs on the contacts the matrices miss
  (care homes, grandparenting), not on susceptibility. The one real age effect on susceptibility is
  subtype: people imprinted on H1 or H2 in childhood, i.e. born before 1968, are at higher risk from
  H3N2 -- a candidate driver for the learning layer rather than a fixed parameter.
- **Transmissibility between countries.** The contact matrices are rescaled away; keeping their
  spectral radius (appendix D) would let demography and mixing set part of each country's R0 and leave
  `S0_c` a purer residual. One line of C++; a model comparison, not a default.
- **Age-specific season visibility.** `delta_s` scales every age group alike, so a season visible
  mainly in the elderly is represented only in aggregate.
- **Spatial spread with outside information.** The dispersion of regional peak times within each
  country would let a fixed per-country spread be set from data rather than fitted.

**About the next step, the learning layer.**
- The season effects on `S0` and on visibility are the quantities to relate to drivers (subtype,
  vaccine match, prior-season burden). On simulated data the fit-then-regress procedure returns a
  known driver effect on both without bias (section 6), so it can be used as it stands.
- With only two seasons standing apart in `S0`, the driver analysis will have little variance in
  susceptibility to explain and a lot in visibility.
- Country-level conclusions should wait for the caveat in point 1 above to be tested -- the
  per-country R0 model comparison is the first candidate.

## 8. The figures, and what to look at in each

All in `output/joint_model/`, regenerated by `run_joint_model.R` (01-15) and `run_joint_recovery.R` (16).

| figure | shows | look at |
|---|---|---|
| 01 data panel | which country-seasons exist, their sources and grid lengths | the gaps, and how short 2025/2026 is |
| 02 data features | the four facts that force the design | panel 1 (levels are countries), panel 2 (seasons co-move) |
| 03 mechanism | what each parameter does to one wave, and what one wave can tell apart | the bottom row: the ties of section 5; the children's modifier panel |
| 04 design | what varies where, and what is fixed | the fixed rows: R0 and the spread; the one global children's modifier |
| 05 arrival | the fitted seeds as arrival times | countries that are consistently early or late |
| 06 fit overview | every country-season, model against data | where the curve misses the peak |
| 07 per country | the three age groups, season by season | the age pattern of the misses |
| 08 noise budget | noise the fit needed against the data's own | the 2.6x gap, per country |
| 09 season `S0` | how easily each season spread | 2015/2016 and 2025/2026 against the flat middle |
| 10 season visibility | how visible each season was per infection | the 2017/2018 peak and the 2023/2024, 2025/2026 troughs |
| 11 country `S0` | how easily influenza spreads in each country | the ranking, and its caveat |
| 12 country reporting | how much of each country's infection is seen | the twentyfold spread; the level is the pin's |
| 13 age reporting | how readily children and the elderly are seen | children about 1.2 times adults, the elderly about 1.8 times |
| 14 attack rates | who gets infected | children highest, the elderly lowest; the pattern, not the level |
| 15 data or prior | contraction per parameter family | nothing near zero |
| 16 recovery | a known truth, simulated and refitted | points on the diagonal; the children's modifier above it, pulled by its prior |

Figures 17 (sensing with `S0` or R0) and 18 (what makes a wave fat) belong to the investigations in
appendices G and E and were produced at commits `7b04681` and `c9286b7`, under the former pin of 1.5.

---

# Appendices: the evidence trail

Each appendix is an investigation as it was run, dated, with the numbers measured at the time. Where
the model has changed since, a note at the top says what applies now.

## Appendix A. What was cut relative to the compartmental pilot, and why

| Cut | Why |
|---|---|
| **The Kalman filter, and with it the process noise `q` and the initial-state covariance `P0`** | Measured: the filter only behaved when `q` was fixed small and `P0` set to zero, and with those the EKF optimum agreed with the deterministic one in every parameter (S0 0.830 vs 0.834, `R0_s` within 0.02). It was costing 7x the runtime and two unfittable knobs to reproduce the answer the deterministic fit already gave. Retained as an option to add back once the innovation diagnostics exist. |
| **The Gaussian observation with count-like variance** | Only needed because a Kalman update requires a Gaussian. Dropping the filter frees the natural choice. At the fitted dispersion the Gaussian put 15% of its predictive mass below zero at any count level, and 29% of observed cells are exactly zero. |
| **The two-stage fit, the engine switch, the age-susceptibility switch, the per-country `R0_s` diagnostic stage** | Machinery that existed to manage the filter or to answer questions now answered (decisions.md, 2026-09-11). |
| **Per-country-season season deviation** | Now one shared value per season (owner decision), which is both the scientific hypothesis and 88 fewer parameters. |
| **Age-specific initial immunity, age-specific susceptibility of the young** | Parked by decision. The young's absolute attack rate is therefore probably too low; recorded as a limitation, not fitted. |

## Appendix B. What was fixed relative to the pilot

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

## Appendix C. What the data layer does and does not support (audited 2026-09-16)

The model code had been audited; the data layer under it had not. It reconstructs correctly — all 86
country-season observation matrices rebuild **exactly** from the raw streams by independent code (max
count difference 0), the age bands and populations are mutually consistent to the person, missing cells
are genuinely missing, and the season boundaries are the documented 1 August. Five things nonetheless
qualify what can be read off the output (a sixth, the positivity-encoding
exclusion, is in the box at the top of this file and in
`documentation/to_confirm_with_surveillance.md`):

- **Season support is uneven, and the thinnest season carries the highest `R0`.** `2025/2026` rests on
  **7 of the 12 countries** (DK, EE, FR, BE, IE, PL, HR), on grids of 34–41 weeks against 52–53
  elsewhere, and 774 observed cells against 1497 for 2023/2024 — and it carries the largest season
  effect of the eight. `jm_summary_season` reports `n_country`, `obs_cells` and the grid range beside
  every season effect, so no season-level number is read without its sample size.
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
  were never offered to the model. Since the season effects are shared across the countries in the
  design, which countries are in it is a substantive choice.

## Appendix D. R0 in a country: what the contact matrix, the age susceptibility and the pyramid do, and do not do

> Written with the elderly susceptibility factor `sigma_eld` in the model. It was removed on 2026-09-27
> (appendix L): the matrix is now used as rescaled, with no age re-weighting, and everything below
> holds with `sigma_eld` = 1.

**R0 is fixed in every country and every season** -- at 2.0 since 2026-09-26, at 1.5 before (owner decisions). This section is
the honest account of what that means once age mixing enters, because it is easy to believe the
model has a country-specific R0 when it does not.

**How the dynamics are built.** Each country has its own contact matrix `C_c` (Prem et al., collapsed
to three groups with that country's population weights, so the pyramid already shapes it). It is
rescaled to spectral radius one. The elderly susceptibility `sigma_eld` then multiplies the elderly
row, and the result is **rescaled to spectral radius one again**. The force of infection uses
`beta = R0 x gamma` on that matrix. Consequence: the realised reproduction number at full
susceptibility is **exactly the pinned R0 in every country**, whatever its matrix, its pyramid, or `sigma_eld`.

**What the matrix, `sigma_eld` and the pyramid therefore control.** *Who* gets infected, not *how
many*. They set the age distribution of infections and the age-specific attack rates — the elderly
attack rate is 17% against 26% for adults (at R0 = 2.0) precisely because of this structure — but they cannot move
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
judgement after the data could not separate the two (appendix G): `S0` is a blunt sensor of
susceptibility and infectivity together, and no parameter in the model carries R0. **It also carries a
country's spatial structure**: a country whose cities are hit at different times has a fatter national
wave than one well-mixed epidemic, and with the spread fixed at none that fatness reads as a lower
`S0` (a slower, broader local wave). A tested spread parameter could not be told apart from `S0` in
practice (appendix E), so this too is by design.

**On the visibility side the age structure IS country-specific, correctly.** Age-specific reporting
(`off_c,young`, `off_c,eld`) varies by country, and the pyramid enters the expected count through the
group populations `N_a`, so a country with more elderly sees more elderly consultations. What the
model does **not** have is age-specific *season* visibility: `delta_s` scales every age group by the
same factor, so an H3N2 season that is disproportionately visible in the elderly is represented only
in aggregate. A natural extension, not a defect.

**The alternative, for a later model comparison.** Skip the second renormalisation. Then the realised
R0 in country `c` is `R0 x rho(D_sigma C_c)`, which varies by country through *known* demography and
mixing — a mechanistically-determined transmissibility — leaving `S0_c` a purer susceptibility
residual. It is a one-line change in the C++ and it changes what `sigma_eld` means (raising it would
then raise overall transmissibility too). Worth fitting side by side once the learning layer starts;
not done now, so that one thing changes at a time.

## Appendix E. Spatial spread, and what makes a wave fat: R0, S0 and tau (2026-09-26)

> Measured at the former pin of R0 = 1.5. Its conclusions carry over to 2.0 unchanged: R0 and `S0` are one
> number to the data (point 1 below), and the pin sweep (appendix F) moves no relative result. Outcome: the
> spread is fixed at none.

**Outcome first (owner, 2026-09-26): the fatness is fixed, not fitted.** A spread parameter clearly
interferes with `S0` (below), so it is pinned like R0 -- at none, `tau_fixed = 0` -- and `S0` is
acknowledged to capture the spatial structure of a country too: the overlay of its local waves, which
makes the national wave fatter than one well-mixed epidemic would be. The spread stays in the code as a
fixed setting, for the day outside information on regional peak timing can set it; the fitted variants
and the analysis below are reproducible from commit `c9286b7` (`run_tau_analysis.R`, figure 18).

**The extension.** A country is not one well-mixed population: its cities are hit at slightly
different times. The model keeps one local epidemic per country-season and treats the country's many
local epidemics as copies of it whose start times are spread normally, sd `tau` days, around the
modelled one. The national incidence is the local incidence convolved with that normal kernel -- on
the daily grid, before the weekly aggregation, each day weighted by the normal mass in its bin and the
kernel cut at 7 sd. As tested, `log tau ~ N(log 7 days, 0.7)`, shared by every country and season or
one per country. Figure 03 has a panel showing what a fixed spread would do.

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

**Decision (owner, 2026-09-26).** The spread is **fixed at none** and `S0` carries the spatial overlay.
A fitted shared spread settles at about two days and changes nothing (`S0` rank correlation 1.00 with
no spread), so fixing it at zero costs nothing measurable. A per-country spread fits better in total
but not wave by wave, where a real country difference would show; its values do not follow country size
and in Estonia fit noise; and within each country it trades against `S0` -- in the simulation where it
is the true model, the country ranking of `S0` came back at 0.69, against 0.94-0.99 without it. **If a
per-country spread is ever to be modelled, it needs outside information to set it**: the observed
dispersion of regional peak times within each country, from sub-national surveillance or published
estimates, as a fixed `tau_fixed` or as a prior. That would take it off the ridge instead of letting it
float there, and it is what more of the same national data cannot do.

## Appendix F. The R0 pin: what its value does (2026-09-26)

> **Decision (owner, 2026-09-26): R0 = 2.0**, with the `S0` prior placed on R0 x `S0` so it moves with the pin.
> Measured before the children's modifier (appendix L). That modifier acts on the logit scale, so the
> children-to-adult ratio it produces depends on the level of `S0` and with it on the pin: its prior
> centre is calibrated at R0 = 2.0, and `run_kappa_profile.R` recalibrates it for any other pin.

R0 and `S0` are exactly one number to the data (appendix E, point 1), so the pin's value
can matter through three things only: the ceiling `S0 <= 1`, the logit-additive form of the season and
country effects, and the ABSOLUTE scale -- at a fixed R0 x S0 the infections scale by 1/R0 and the
reporting fraction by R0. `run_pin_sweep.R` refits the working model at five pins, the `S0` prior
centre moving with the pin (0.75 x 1.5/R0) so the prior sits on the same R0 x S0 everywhere:

| R0 pin | log-lik vs 1.5 | `S0` median / max | `S0` > 0.95 | R0 x `S0` median | country `S0` rank vs 1.5 | attack rate, all ages / adults / 65+ | reporting, adults |
|---|---|---|---|---|---|---|---|
| 1.5 | 0 | 0.86 / 0.985 | 8 | 1.29 | 1 | 32% / 35% / 22% | 3.1% |
| 1.7 | +9.5 | 0.76 / 0.90 | 0 | 1.29 | 0.99 | 28% / 31% / 20% | 3.5% |
| 2.0 | +11.6 | 0.64 / 0.78 | 0 | 1.29 | 0.99 | 24% / 26% / 17% | 4.1% |
| 2.5 | +11.6 | 0.51 / 0.63 | 0 | 1.28 | 0.99 | 19% / 21% / 13% | 5.1% |
| 3.0 | +11.2 | 0.43 / 0.53 | 0 | 1.28 | 0.99 | 16% / 17% / 11% | 6.2% |

- **Off the ceiling the pin changes nothing relative.** From 1.7 up the fit is flat within 2 nats, the
  country ranking of `S0` holds at 0.99, the season effects correlate at 0.996 or more between pins and
  season visibility at 0.998. None of the differences is evidence with the wave as the unit (sign test
  p 0.39-1.0, every bootstrap interval spans zero): the gain from 1.5 to 1.7 is the ceiling releasing
  Croatia and 2025/2026, not a signal about transmissibility.
- **1.5 is the only pin that differs, through the ceiling alone**: 8 cells press 1, and 2025/2026 is
  compressed (a typical country's R0 x `S0` 1.38 against 1.40 at every higher pin). The six seasons
  between 2015/2016 (lowest, 1.23) and 2025/2026 (highest) lie within 1.27-1.30 of each other, so their
  RANKS swap between pins; their values do not move.
- **What the pin does set is the absolute scale**: attack rates fall as 1/R0 and reporting fractions
  rise as R0. R0 x `S0` itself stays at 1.29 -- which is where a systematic review puts the
  reproduction number of seasonal influenza (median 1.28; Biggerstaff et al. 2014), so the data and the
  literature agree on the product. How it splits into R0 and `S0` is the pin's choice; for populations
  with little prior immunity the same review's pandemic estimates run from about 1.5 to 1.8.

The recommendation to the owner and the decision are in `documentation/decisions.md` (2026-09-26).

## Appendix G. Sensing with S0 or with R0: the model comparison (2026-09-25)

> Measured at the former pin of R0 = 1.5; the variants it compares are reproducible from commit `7b04681`.

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

## Appendix H. Identifiability of the fixed-R0 model at R0 = 1.5 (2026-09-25)

> Superseded for the current numbers by section 5, which repeats every measurement on the R0 = 2.0 fit; kept for
> the misspecification arms as first run on this model.

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

## Appendix I. Implementation: how it is fitted, and the two guards on the optimiser

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

## Appendix J. Recovery under the previous model (measured 2026-09-12)

> Measured under the model with a fitted per-season R0 on the 86-cell design. Section 6 has the recovery of
> the current model.

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

## Appendix K. Are the headline claims robust? (stress-tested 2026-09-14)

> Measured under the model with a fitted per-season R0. The first result's analogue for the current model --
> seasons differ far more in visibility than in how easily they spread -- is in section 4, and the noise-gap
> result stands as measured.

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

## Appendix L. Age and susceptibility: the evidence, the children's modifier, and the elderly factor removed (2026-09-27)

**The owner's questions.** Would expert knowledge say that children and the elderly both start a season
more susceptible than adults? And do cohort studies show that an older immune system raises the risk of
infection, or only of severe outcomes?

**What the literature says.**

- **Children are infected most, because they carry the least immunity.** Cohorts that test everyone,
  whatever their symptoms, find children infected most often: PHIRST (South Africa 2016-2018,
  twice-weekly PCR; Cohen et al., Lancet Glob Health 2021) about 1.7 times as often as adults once its
  age bands are collapsed to ours (1.65-1.80; decisions.md, 2026-09), and the Tecumseh study (Monto et
  al., Am J Epidemiol 1985) highest in children and lowest in older adults. Household studies find
  children more susceptible to infection (Tsang et al., Trends Microbiol 2016). The mechanism is prior
  immunity: antibodies to influenza accumulate over a lifetime of infections (Kucharski et al., PLoS
  Biol 2015), so children enter each season with the least.
- **Older immune systems raise severity, not infection.** About 90% of influenza-associated deaths are
  in people aged 65 and over (Thompson et al., JAMA 2003), and their odds of an antibody response to
  vaccination are a quarter to a half of young adults' (Goodwin et al., Vaccine 2006). But the same
  cohorts find them infected no more often than adults: 0.80 times in PHIRST, least of all ages in
  Tecumseh. Decades of exposure leave cross-reactive immunity: in 2009 about a third of people over 60
  already carried antibodies to the new pandemic H1N1, and almost no children did (Hancock et al.,
  NEJM 2009).
- **The exception is subtype, through imprinting.** The first influenza met in childhood shapes lifelong
  protection by haemagglutinin group (Gostic et al., Science 2016): people born before 1968, imprinted on
  H1 or H2, are at higher risk from H3N2 and lower from H1N1 (Gostic et al., PLoS Pathog 2019). That is
  an age effect whose sign depends on the season's subtype -- a driver for the learning layer, not a
  constant.

**What that means for the model.**

- *Children.* The model had no way to express it: `S0` was the same for every age, and a per-contact
  susceptibility for children had been rejected by the wave shapes (decisions.md, 2026-09, which named
  children's lower initial immunity as the candidate left to test). Children were infected 0.93 times as
  often as adults, and their whole excess of ILI+ went into a reporting offset of 1.95.
- *The elderly.* `sigma_eld` (2.2 per contact) could not be higher susceptibility to infection; at most
  it was exposure the contact matrices miss. The owner removed it (below).

**The modifier.** `logit S0_young = logit S0_{c,s} + kappa`, `kappa = exp(log kappa) > 0`: one number for
all countries and seasons (biology), able only to raise children's `S0` (the owner's constraint), and
never past 1 (the logit). One parameter. As `kappa` goes to zero it is the previous model:
against the engine of commit `da3a80e` the working fit's log-likelihood agrees to 6.5e-9 nats and every
fitted mean to 7e-13, and the C++ agrees with the base-R reference to 1e-10.

**What the data say** (`run_kappa_profile.R` on the current model, without the elderly factor: `kappa`
held at each value, everything else refitted from two starts -- the working fit and the nearest value
already profiled -- because from one start the fit at zero once stuck 45 nats worse, Italy 2024/2025's
wave moved early and its ERVISS baseline absorbing the counts. Medians over country-seasons. With the
elderly factor in the model the profile was the same within 2.5 nats.)

| `kappa` | children's / adults' `S0` | log-lik vs none | children-to-adult infection | children's reporting x adult |
|---|---|---|---|---|
| 0 | 0.64 / 0.64 | 0 | 0.94 | 1.87 |
| 0.25 | 0.68 / 0.63 | -2.4 | 1.08 | 1.64 |
| 0.5 | 0.73 / 0.62 | -5.6 | 1.20 | 1.48 |
| 1 | 0.80 / 0.59 | -15.4 | 1.45 | 1.22 |
| 1.5 | 0.85 / 0.57 | -21.8 | 1.70 | 1.04 |
| 2 | 0.90 / 0.54 | -27.9 | 1.95 | 0.92 |
| 3 | 0.96 / 0.52 | -37.6 | 2.38 | 0.77 |

- **Surveillance alone prefers no modifier, weakly.** As `kappa` rises the children's reporting offset
  falls in step: the level of children's counts is held and only its attribution moves (the tie of
  section 5). What could break the tie is shape -- more susceptible children drive the epidemic more, so
  their wave should lead and steepen relative to the adults' -- and the data barely show it, the verdict
  the per-contact test also gave (decisions.md, 2026-09).
- **With the wave as the unit, the cost is not evidence at the values that matter.** At `kappa` = 1, 36
  of 85 waves favour the modifier (sign test p 0.19, Wilcoxon 0.14) and the bootstrap intervals over
  waves and countries span zero; at 1.5 likewise (37 of 85; p 0.28 and 0.10). Small values carry a small
  consistent cost (at 0.5, -5.6 nats, Wilcoxon p 0.02: the wave bootstrap excludes zero, the country
  bootstrap does not), and so do large ones: `kappa` = 2 and 3, children infected 2-2.4 times adults,
  are rejected at the margin (Wilcoxon p 0.03, sign p 0.05).
- **The cost is two countries'.** At `kappa` = 1, Ireland carries -12.0 and Croatia -5.6 of the -15.4;
  the other ten together favour the modifier slightly (+2.2), Italy, France, Denmark and Poland most.

**The decision (owner's instruction, 2026-09-27): the cohort evidence anchors it.** The data cannot
separate `kappa` from children's reporting and do not reject the cohort value, so the prior carries it:
`log kappa ~ N(log 1.5, 0.2)`. The centre is where the refitted model reproduces PHIRST's ratio of 1.7
(1.49 by interpolation in the profile; 1.59, and a centre of 1.6, while the elderly factor was in the
model); the sd admits ratios of about 1.45-2.0, for the uncertainty in carrying a South African cohort
to Europe. The fit: `kappa` = 0.94 (95% 0.70-1.28), children infected 1.42 times as often as adults
(interquartile over country-seasons 1.30-1.60), contraction 0.23 -- the data pull it partway back and
the prior does most of the work, by design. Against the same model with `kappa` held at zero it costs
14.6 raw nats.

**What moved with it, and what did not** (measured when the modifier was added, with the elderly factor
still in the model).

| quantity | without the modifier | with it |
|---|---|---|
| attack rate, children / adults / elderly (median) | 25% / 26% / 17% | 35% / 23% / 16% |
| children's reporting x adult (median over countries) | 1.95 | 1.25 |
| reproduction number at season start (median) | 1.29 | 1.29 |
| season effects on `S0`; visibility | | Pearson 0.999; 0.999 with the fit without it |
| country `S0` ranking | | Spearman 0.90 with the fit without it; Ireland moves to last |
| country `S0` against data quality (Spearman with the dispersion) | 0.75 | 0.52 |

**Caveats.**

- PHIRST is a South African cohort: a younger population, higher HIV prevalence, larger households.
  Flu Watch (England; Hayward et al., Lancet Respir Med 2014), which measured infection by serology at
  all ages, is the natural European anchor if its age-specific rates can be extracted.
- `kappa` acts on the logit scale, so the ratio a given value produces depends on the level of `S0`,
  and so on the R0 pin and the age structure: the prior centre is calibrated at R0 = 2.0 without an
  elderly factor and must be recalibrated when either changes (`run_kappa_profile.R` prints the new
  centre).
- In the recovery study (section 6) the modifier comes back at its prior's pull, not at the truth --
  simulated at 0.94, estimated at 1.22 -- and the level of adults' `S0` shifts with it. That is the
  design working as intended: the modifier, and with it the level of `S0`, is the cohort's; the ranking
  and the season effects are the data's.

**The elderly factor removed (owner, 2026-09-27).** On the evidence above and on parsimony, `sigma_eld`
is gone: its slot, its prior and the re-weighting of the contact matrix, so 181 parameters remain and
each matrix is used as rescaled. With the factor at 1 the new engine reproduces the old to 1.2e-8 nats,
and the refit reaches the same optimum from a cold start and from the old fit (to 0.000 nats).

| | with the elderly factor | without |
|---|---|---|
| log-likelihood (parameters) | -50800.3 (182) | -50811.1 (181) |
| waves fitted better | 58 of 85 | 27 of 85 |
| elderly attack rate, median (x adults) | 16% (0.69) | 6.5% (0.29); PHIRST 0.80 |
| elderly reporting x adult, median over countries | 0.78 | 1.83 |
| children's modifier `kappa` | 1.02 (prior centre 1.6) | 0.94 (prior centre 1.5) |
| season effects, visibility, country `S0` ranking | | correlation 1.000 with the fit with the factor |

- **The data would keep it, modestly but consistently.** 10.9 raw nats is small -- 0.13 per wave -- but
  it is evidence with the wave as the unit: 58 of 85 waves fit better with the factor (sign test p 0.001,
  Wilcoxon p 0.0001), every bootstrap interval excludes zero (waves -17 to -5 nats, countries -16 to -6,
  seasons -19 to -3), and 11 of 12 countries lean the same way. Half the cost is in the elderly's own
  counts and half in the children's: the elderly's share of transmission shapes the whole wave a little.
- **What it did.** It moved the elderly's excess of cases between infection and visibility. With it,
  the elderly were infected 0.69 times as often as adults and seen 0.78 times as readily per infection;
  without it, 0.29 times and 1.8 times. The first matches the cohort's infection ratio; the second
  matches the literature's reading, since more severe illness makes an infection likelier to reach a
  doctor. Neither touches what the model is for.
- **Reading.** The literature rules out higher susceptibility, not higher exposure: the contact matrices
  may under-record the elderly's contacts (care homes, grandparenting), and PHIRST's elderly live in
  multigenerational households, so how far below adults the European elderly are infected is not known
  well enough to fix a correction. The removal trades a small, consistent loss of fit and a lower
  elderly attack rate for one parameter fewer and a model with nothing to reinterpret; nothing the
  learning layer uses changes.
