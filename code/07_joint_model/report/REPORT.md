::: {.titleblock}
[Status report · 27 September 2026]{.kicker}

# Seasonal influenza in twelve EU/EEA countries {.title}

[A joint transmission model of eight seasons: the data, the evidence, the model, its behaviour and its fit]{.subtitle}

[Project flu_seasonal_variation_analysis · repository reneniehus/flu-susceptible-fit, branch flu_comp_model · model code in code/07_joint_model/]{.meta}
:::

::: {.keynums}
::: {.kn}
[12 × 8]{.n}
[countries × seasons; 85 country-seasons observed]{.l}
:::
::: {.kn}
[8,685]{.n}
[weekly counts in three age groups, 31% of them zero]{.l}
:::
::: {.kn}
[181]{.n}
[parameters: 15 shared across Europe, 166 local to one country]{.l}
:::
::: {.kn}
[0.91]{.n}
[median correlation of fitted and observed weekly counts]{.l}
:::
:::

::: {.abstract}
**Abstract.** We fit weekly influenza-positive ILI consultations (ILI+) from twelve EU/EEA countries over eight seasons (2014/15–2018/19 and 2023/24–2025/26) with one age- and vaccination-structured SIR model per country-season and a negative-binomial observation model, estimated jointly in 181 parameters. The model is best read as a structured sensor of wave shape. It separates how easily influenza spread — a susceptibility-like index S~0~ with one effect per country and one per season shared across Europe — from how visible that spread was: a reporting level per country and a visibility per season.

Three inputs are set from outside the data, because the data cannot separate them. R~0~ is pinned at 2.0, since the data fix only R~0~ × S~0~. Spatial spread is fixed at none. Children start each season more susceptible, by an amount anchored on cohort evidence (PHIRST). The elderly carry no susceptibility factor: ageing raises the severity of influenza, not the risk of infection.

Seasons differ little in how easily they spread — a typical country's adult R~0~ × S~0~ ranged from 1.12 to 1.30 — and much more in visibility, from 0.57 to 1.84 times normal. The two move against each other (correlation −0.61). The share of adult infections that is reported ranges from 1.0% to 24.5% between countries. Children are infected 1.4 times as often as adults.

The model reproduces the waves: the median correlation is 0.91, and the peak falls within one week in 58 of 85 country-seasons. Its limits:

- It needs 2.6 times more noise than the data's own week-to-week scatter.
- It cannot follow a country whose season departs from the European pattern.
- It ranks countries reliably only if transmissibility is equal across them.

On simulated data it recovers the season effects, visibility and the country ranking (0.97–0.99), and a driver's effect on the season parameters without bias. The season effects are ready for the learning layer.
:::

## 1. Why a joint model

The project asks why influenza seasons differ. Raw ILI+ curves cannot answer that directly. A large season may be one in which more people were infected, or one in which more of the infections were seen — and every surveillance system sees a different share. The joint model separates the two for every season and every country at once, borrowing strength across Europe. It hands the season effects to the learning layer, which relates them to drivers such as the dominant subtype, vaccine match and the burden of the previous season.

## 2. The data

The data are weekly ILI+ per 100,000 in three age groups (0–14, 15–64 and 65+) from ECDC surveillance, stitched week by week into one panel. RespiCompass supplies the seasons up to 2023/24. ERVISS supplies them from 2024/25, as ILI rate × sentinel positivity, or non-sentinel positivity for Croatia. Seasons run from 1 August, and the COVID seasons 2019/20–2022/23 are omitted.

Rates are converted to counts with each age group's population, a scale device the reporting parameters absorb. A country enters with at least five age-complete seasons, and twelve do (Figure 1). Iceland, Malta and Austria would also qualify and are not included. Because the season effects are shared, which countries are in the design matters.

One country-season is excluded provisionally. ERVISS encodes a week with tests but no detections row ambiguously: it may be a zero or an unpublished count. Fourteen weeks of Czechia 2024/25 cannot plausibly be zero. The question is open with the surveillance team.

![**Figure 1. The panel.** One tile per country and season; the number is how many weeks carry an observation, the colour the source. Blank tiles: no usable data. The cross marks Czechia 2024/25, excluded pending the ERVISS encoding question.](fig01_data_panel.png)

Four features of the data dictate the design (Figure 2):

1. **The level is a country property.** Peak ILI+ differs about 150-fold between countries, and each country keeps its rank across seasons. That reflects how much a surveillance system sees, not how many people were infected. So each country has its own reporting level.
2. **Seasons move together.** Once each country's level is divided out, a big season is big almost everywhere. So each season has one effect on S~0~ and one on visibility, shared by all countries.
3. **About a third of the counts are zero.** So observations are negative-binomial counts, not Gaussian.
4. **Waves arrive weeks apart:** peaks span 14 weeks. So every country-season has its own seed, the model's only handle on timing.

![**Figure 2. Four features of the data, and the model choice each forces.** (1) Peak weekly ILI+, one point per season (log scale). (2) Each country's peak relative to its own median; grey lines are countries, orange the median. (3) Weekly counts in one age group. (4) Week of the observed peak, all 85 country-seasons.](fig02_data_features.png)

## 3. The model

### How it works

Each country-season is an independent epidemic. It is an SIR model in three age groups, each split into vaccinated and unvaccinated, integrated daily from a seed on 1 August to a 53-week horizon. The 65+ group is vaccinated in one pulse on 1 October at the reported national coverage. Vaccination reduces infection by 25%, ILI given infection by 20% and onward transmission by 20%. The force of infection on age group *a* is

::: {.eq}
λ~a~(t) = β Σ~j~ C~aj~ [ I~j~^unvacc^(t) + 0.8 I~j~^vacc^(t) ],  β = R~0~ γ.
:::

Here C is the country's contact matrix (Prem et al. [12], built from POLYMOD surveys [11]), rescaled to spectral radius one. This makes the reproduction number at full susceptibility exactly the pinned R~0~ = 2 in every country. The infectious period 1/γ is 3.6 days.

Infection per contact is the same at every age. What differs is the susceptible fraction at the start of the season: S~0~ for adults and the elderly, and a higher value for children:

::: {.eq}
logit S~0,c,s~ = S~0,c~ + x~s~   and   logit S~0,c,s~^children^ = logit S~0,c,s~ + κ,  κ > 0.
:::

Each week's new infections in age group *a* become counted consultations through the reporting chain

::: {.eq}
E[y~a~] = new infections~a~ × c~c~ × δ~s~ × o~c,a~ + b~c,source~,  y~a~ \~ negative binomial with dispersion φ~c~.
:::

| Parameter | Varies by | Count | Meaning |
|:----------------------|:------------------------------------------|:----|:------------------------------------------------------------------|
| S~0,c~ | country | 12 | country's susceptible fraction at the average season (logit scale) |
| x~s~ | season, shared by all countries; averages zero | 7 | season effect on logit S~0~ |
| κ | nothing: one number, only positive | 1 | how much more susceptible children start a season |
| c~c~ | country | 12 | share of adult infections that become ILI+ consultations |
| δ~s~ | season, shared; averages one | 7 | season visibility: consultations per infection relative to normal |
| o~c,a~ | country × age | 24 | how readily children and the elderly are seen, relative to adults |
| φ~c~ | country | 12 | negative-binomial dispersion |
| b~c,source~ | country × data source | 21 | small off-season floor |
| I~0,c,s~ | country-season | 85 | the seed on 1 August, which sets when the wave arrives |
| fixed | — | — | R~0~ = 2.0; spatial spread none; 1/γ = 3.6 days; vaccine effects; contact matrices |

: **Table 1. The parameters.** 181 in all; 15 are shared across countries (7 + 7 season effects and κ).

The likelihood is separable given the shared parameters, so it is maximised by block coordinate descent. All countries' local blocks are optimised in parallel, then the shared block, repeatedly, followed by a joint polish. Local blocks are multi-started, and a guard stops any country collapsing to a flat line.

Intervals are 95% intervals from the curvature of the log-posterior. The log-posterior is written in C++ and checked against an independent base-R implementation to 10^−10^.

![**Figure 3. What varies where.** Each row is a parameter; the squares say whether it may differ between seasons, countries or age groups, and the right column how many numbers it contributes.](fig05_design.png)

### How to read the parameters

- **S~0,c,s~**: how easily influenza spread in that country and season, expressed as a susceptible fraction at R~0~ = 2. It is a deliberately blunt sensor. Susceptibility, strain transmissibility, demography and mixing (the contact matrix's scale is removed) and the overlay of a country's local waves all land in it. Read its ranking; its level is set by the R~0~ pin.
- **x~s~**: this season, for all of Europe.
- **c~c~ and δ~s~**: how many ILI+ consultations one infection produces in that country, and how many more or fewer in that season. That covers how symptomatic the strain was, how many consulted and how many were tested. The data see only their product, so the constraint that δ averages one separates them.
- **κ**: how much less immunity children carry into a season than adults. This is biology, one number for Europe.
- **o~c,a~**: surveillance — how readily each age group is seen, per country.
- **I~0,c,s~**: timing only (Figure 14).
- **φ~c~**: the country's noise. It also absorbs misfit, which makes it the adequacy diagnostic (Section 6).

### How the model behaves

Figure 4 moves one parameter at a time on one fitted wave.

- **Susceptibility (season and country effects on S~0~):** the two are the same lever. A higher S~0~ makes the wave rise faster, peak earlier and grow taller.
- **Reporting and visibility:** these scale the curve without changing its shape.
- **The seed:** moves the wave sideways and changes nothing else.
- **κ:** raises and advances the children's wave.
- **Spatial spread:** would lower and widen the peak; it is fixed at none.

![**Figure 4. What each parameter does to a wave.** One fitted country-season (Poland 2024/25), adults unless stated. Dark grey: the fit. Blue and orange: one parameter moved down and up, all else held. The settings are: effects on S~0~ ±0.25 logit; seed ÷100 and ×100; reporting and visibility −40% and +60%; 65+ reporting halved and doubled; κ none and doubled; spread none and 14 days. Every curve is produced by the model's own code.](fig03_mechanism.png)

What one wave can and cannot tell apart decides the design (Figure 5):

- **More susceptible or more visible:** a bigger season from more susceptibility rises faster and peaks earlier than one from more visibility. The shapes differ, so the data separate them.
- **Less susceptible or more spread out:** a fatter wave from lower S~0~ and one from spatial spread differ only in the low tails. So the spread is fixed rather than fitted.
- **Reporting level and season visibility:** these are exactly one number. Only the average-one constraint on visibility separates them.
- **R~0~ and S~0~:** these are exactly one number too. Any factor on R~0~ is undone by S~0~, the seed and reporting, hence the pin.

![**Figure 5. What one wave can and cannot tell apart** (log scale, so exponential growth is a straight line). Left: a bigger season from more susceptibility (blue) or more visibility (orange), tuned to the same peak. Middle: a fatter wave from lower S~0~ (green) or from spatial spread (violet), tuned to the same peak and width. Right: reporting up 50% with visibility down by the same factor gives the identical curve.](fig04_ties.png)

## 4. What informs the parameters: studies and data

| Quantity | Set by | Evidence |
|:----------------------|:------------------------|:-------------------------------------------------------------------|
| R~0~ = 2.0 | judgement (a pin) | The data fix only R~0~ × S~0~: the reproduction number at season start is 1.28 (median over country-seasons), where a systematic review puts seasonal influenza (median 1.28) [1]. The pin sets the absolute scale — attack rates scale as 1/R~0~, reporting as R~0~ — and no relative result; 2.0 keeps S~0~ well below 1. |
| Spatial spread: none | judgement | A fitted spread traded against S~0~ and cost the country ranking its identifiability; S~0~ absorbs the overlay of a country's local waves. |
| 1/γ = 3.6 days | project notes | Consistent with the serial interval of influenza, about 3.6 days [3]. |
| Vaccine effects | project notes | 0.25 against infection, 0.20 against ILI given infection, 0.20 against onward transmission; 65+ coverage as reported; no coverage below 65. |
| Contact matrices | Prem et al. [12] | Synthetic matrices from POLYMOD [11], collapsed to three age groups; Norway uses the EU average. Only their structure is used. |
| κ | cohort prior, log κ \~ N(log 1.5, 0.2) | PHIRST [2], Tecumseh [10], household studies [14], antibody life course [9]. The data alone prefer less (Figure 6). |
| No elderly factor | literature and parsimony | Severity, not infection [2, 6, 7, 10, 13]. The data would keep a modest correction (Table 3). |
| All other parameters | the data | Weak priors; the data determine them (contraction 0.78–0.97; Figure 7). |

: **Table 2. Where each input comes from.**

### The age evidence

- **Children are infected most, because they carry the least immunity.** Cohorts that test everyone, whatever their symptoms, find children infected most often:
  - about 1.7 times as often as adults in PHIRST [2], once its age bands are collapsed to ours;
  - most often of all ages in Tecumseh [10].

  Household studies find children more susceptible [14], and antibodies accumulate over a lifetime of infections [9].
- **An older immune system raises severity, not infection.**
  - About 90% of influenza-associated deaths are in people aged 65 and over [13], and the elderly respond less well to vaccination [6].
  - Yet the same cohorts find them infected no more often than adults: 0.80 times in PHIRST, and least of all ages in Tecumseh.
  - Decades of exposure leave cross-reactive immunity: in 2009 about a third of people over 60 already carried antibodies to the new pandemic virus [7].
- **The exception is subtype.** Childhood imprinting sets lifelong protection by haemagglutinin group, so people born before 1968 are at higher risk from H3N2 [4, 5]. The sign of this effect depends on the season's subtype, so it is left to the learning layer.

### Children: the surveillance data are nearly silent, so the cohort decides

Surveillance counts alone cannot tell children's susceptibility from how readily children are seen. Raising κ lowers the children's reporting offset in step. The shape signal that could separate the two — children's waves leading and driving the adults' — is weak. The profile prefers no modifier (Figure 6a), but the preference is not significant when each wave counts as one observation:

- at κ = 1, 36 of 85 waves favour the modifier (sign test p = 0.19);
- every bootstrap interval spans zero;
- the cost falls almost entirely on Ireland and Croatia.

So the prior carries κ: log κ \~ N(log 1.5, 0.2). It is centred where the refitted model reproduces PHIRST's ratio of 1.7 (Figure 6b). Its width admits ratios of about 1.45–2.0, allowing for the uncertainty of carrying a South African cohort over to Europe. The fit pulls κ partway back to 0.94 (0.70–1.28), at which children are infected 1.42 times as often as adults.

![**Figure 6. What informs the children's modifier κ.** Each blue point is a full refit with κ held at that value. (a) Log-likelihood relative to no modifier: the surveillance data prefer none, weakly. (b) The refitted model's children-to-adult attack-rate ratio (median over country-seasons) against the PHIRST cohort's 1.65–1.80. Grey band: the prior's 95% range; orange: the fit and its 95% interval.](fig06_kappa.png)

### The elderly: no factor

An earlier version let the elderly be infected 2.2 times as often per contact as the contact matrix implies. The literature rules that out as susceptibility, and the factor was removed for parsimony. The data would have kept it, modestly: the fit is 10.9 nats better with it, and 58 of 85 waves and 11 of 12 countries favour it (sign test p = 0.001; Table 3). The factor only moved the elderly's cases between infection and visibility. Nothing the learning layer uses changed.

| | With the elderly factor | Without (current model) |
|:---------------------------------------------|:--------------------------|:------------------------------|
| Log-likelihood (parameters) | −50,800.3 (182) | −50,811.1 (181) |
| Waves fitted better | 58 of 85 | 27 of 85 |
| Elderly attack rate, median (× adults) | 16% (0.69) | 6.5% (0.29); PHIRST: 0.80 |
| Elderly reporting × adults (median over countries) | 0.78 | 1.83 |
| Season effects, visibility, country ranking | | correlation 1.000 with the other fit |

: **Table 3. The elderly factor, with and without.**

### The data decide the rest

Every other parameter family is determined by the data rather than by its weak prior (Figure 7). Only κ is prior-dominated, by design.

![**Figure 7. The data or the prior?** Prior-to-posterior contraction per parameter family: 0 means the fit returned the prior unchanged, 1 that the data determined it. Dot size: number of parameters.](fig07_data_or_prior.png)

## 5. What the model learned

### Seasons

Seasons differ little in how easily influenza spread and a lot in how visible it was (Figure 8, Table 4).

- **How easily it spread:** only 2015/16 (below) and 2025/26 (above) stand apart in S~0~. The other six lie within the data's resolution of each other.
- **How visible it was:** visibility ranges from 0.57 to 1.84 times normal. On a common log scale that is about ten times the spread of the susceptibility effect (sd 0.40 against 0.04).
- **The two move against each other:** the correlation across seasons is −0.61, and it is not forced by the model (their posterior correlation is 0.23). The seasons that spread most easily turned the fewest infections into consultations.
- **2025/26** carries both extremes on the thinnest data: 7 countries, with series not yet complete.

| Season | Countries | Adult S~0~, typical country | Adult R~0~ × S~0~ | Visibility × normal (95%) | Median attack rate |
|:--------|------:|-----------:|---------:|------------------:|---------:|
| 2014/15 | 11 | 0.588 | 1.18 | 1.40 (1.28–1.54) | 21% |
| 2015/16 | 12 | 0.562 | 1.12 | 1.16 (1.05–1.27) | 16% |
| 2016/17 | 12 | 0.597 | 1.19 | 0.93 (0.85–1.02) | 21% |
| 2017/18 | 12 | 0.585 | 1.17 | 1.84 (1.68–2.02) | 20% |
| 2018/19 | 11 | 0.593 | 1.19 | 1.09 (0.99–1.21) | 22% |
| 2023/24 | 11 | 0.592 | 1.18 | 0.58 (0.52–0.65) | 19% |
| 2024/25 | 9 | 0.598 | 1.20 | 0.99 (0.88–1.11) | 20% |
| 2025/26 | 7 | 0.648 | 1.30 | 0.57 (0.50–0.65) | 29% |

: **Table 4. The seasons.** A typical country is the median country level plus the season effect. Attack rates are all ages, over a full season.

![**Figure 8. The season effects, shared by all countries.** (a) Adult S~0~ of a typical country; shaded: the prior's 95% range, dashed: its centre. (b) Visibility, consultations per infection relative to the average season (log scale). Bars: 95% intervals.](fig08_seasons.png)

### Countries

- **How easily influenza spreads:** Croatia's waves spread most easily (adult R~0~ × S~0~ 1.40), Ireland's and Estonia's least (1.11 and 1.13). The ranking is sharp, with intervals of ±0.01–0.03 on S~0~. It is also conditional: it assumes equal transmissibility and spatial structure across countries (Section 7).
- **How much is reported:** the share of adult infections that becomes an ILI+ consultation runs from 1.0% in Croatia to 24.5% in Italy. The order and spread are the data's; the level is the pin's.

| Country | Adult S~0~ (95%) | R~0~ × S~0~ | Adult reporting | Children × adults | 65+ × adults | Dispersion | Noise excess | Median attack |
|:---|:----------------|---:|---:|---:|---:|---:|---:|---:|
| HR | 0.701 (0.690–0.713) | 1.40 | 1.0% | 2.23 | 1.08 | 1.79 | 2.5× | 36% |
| ES | 0.649 (0.636–0.661) | 1.30 | 2.6% | 1.55 | 1.78 | 1.33 | 2.0× | 27% |
| FR | 0.624 (0.608–0.639) | 1.25 | 7.5% | 1.16 | 1.81 | 0.60 | 2.6× | 25% |
| CZ | 0.623 (0.611–0.635) | 1.25 | 2.1% | 2.24 | 1.46 | 0.37 | 2.5× | 24% |
| DK | 0.605 (0.587–0.623) | 1.21 | 3.1% | 1.33 | 2.78 | 0.40 | 2.4× | 22% |
| BE | 0.600 (0.582–0.617) | 1.20 | 14.2% | 0.75 | 1.85 | 0.29 | 3.9× | 21% |
| IT | 0.592 (0.567–0.616) | 1.18 | 24.5% | 1.07 | 1.90 | 0.98 | 3.5× | 19% |
| NO | 0.587 (0.573–0.600) | 1.17 | 6.1% | 0.24 | 1.28 | 0.65 | 2.3× | 19% |
| PL | 0.583 (0.570–0.596) | 1.17 | 18.3% | 2.73 | 1.93 | 0.15 | 3.0× | 18% |
| NL | 0.574 (0.551–0.597) | 1.15 | 3.4% | 1.08 | 8.63 | 0.18 | 2.5× | 17% |
| EE | 0.565 (0.551–0.578) | 1.13 | 6.8% | 2.77 | 1.25 | 0.24 | 3.4× | 15% |
| IE | 0.557 (0.529–0.586) | 1.11 | 1.4% | 0.34 | 3.40 | 0.75 | 3.0× | 19% |

: **Table 5. The countries**, ordered by S~0~. "Children × adults" and "65+ × adults" are reporting offsets per infection. Noise excess is the fitted noise over the data's own scatter (Figure 13).

![**Figure 9. The countries.** (a) Adult S~0~ at the average season; shaded: the prior's 95% range. (b) Share of adult infections reported as ILI+ (log scale). Bars: 95% intervals.](fig09_countries.png)

### Age

- **Who is infected:** in the median country-season, children 34%, adults 22%, the elderly 6.5%. That is 1.42 times adults for children (PHIRST: 1.7) and 0.29 times for the elderly (PHIRST: 0.80; see Section 6).
- **Who is seen, per infection:**
  - a child 1.24 times as readily as an adult (median over countries; from 0.24 in Norway to 2.8 in Estonia);
  - an elderly person 1.83 times (from 1.1 in Croatia to 8.6 in the Netherlands), consistent with more severe illness bringing the elderly to a doctor.
- **Levels are the pin's:** the median of 20% across all ages would read 27% at R~0~ = 1.5 and 14% at 3.

![**Figure 10. Who is infected, and who is seen.** (a) Modelled attack rate over a full season; thin lines are countries, thick the median of each age group. (b) Age reporting offsets relative to adults, per country (log scale; bars are 95% intervals).](fig10_age.png)

## 6. How well the model describes the data, and where it does not

::: {.full}
![**Figure 11. Every country-season, model against data.** Points: observed ILI+ per 100,000 (age groups pooled by population); lines: the fitted model. One susceptibility effect and one visibility per season must serve every country at once, so a panel that misses carries information about that shared assumption.](fig11_fit_overview.png)
:::

The model reproduces the waves (Figure 11):

- **Shape:** fitted and observed weekly counts correlate at a median of 0.91 per country-season (interquartile 0.85–0.95). By age group the medians are 0.88 for children, 0.91 for adults and 0.88 for the elderly.
- **Timing:** the fitted peak falls within one week of the observed peak in 58 of 85 waves, and within two weeks in 70.
- **Height:** right on average (observed/fitted median 0.92) but scattered. The model is more than 25% too low in 18 waves and more than 25% too high in 29.

Figure 12 shows one country by age. The age profile is broadly right. Where the model misses, it usually misses all ages together (Denmark 2014/15 too high, 2015/16 too low): a season departing from the European pattern rather than a fault in the age structure. The one clear age-specific miss is the elderly in 2023/24.

![**Figure 12. One country by age (Denmark).** Points observed, lines fitted, per age group and season. The age profile comes from the contact matrix, κ and Denmark's two reporting offsets, none of which varies by season.](fig12_fit_DK.png)

**Where the model does not describe the data.**

1. **Countries whose season departs from Europe's.** One season effect serves every country, so the model cannot follow a country-season that breaks the European pattern. Italy 2023/24 and Ireland 2015/16 peaked more than twice as high as the model; Norway 2023/24 and Croatia 2015/16 at about a third. A country-specific season deviation, shrunk towards zero, is the natural extension.
2. **Within-season wander** (Figure 13). The fit needs 2.6 times the data's own week-to-week scatter (median; 2.0–3.9 by country; 2.1–3.7 however the scatter is measured). The deterministic curve cannot follow sharp peaks, shoulders and second humps (for example France 2017/18 and Estonia 2023/24), and writes that misfit off as noise. Candidate remedies are process noise, a second strain and a school-holiday contact effect.
3. **Estonia.** Its series are the noisiest: five of the 15 waves whose fitted peak is three or more weeks off, and four of the five lowest correlations (0.35–0.64).
4. **Timing through the seed** (Figure 14). Every wave is seeded on 1 August, so its seed size sets its timing. Croatia's fast waves need seeds near 10^−11^, three prior standard deviations below the rest. Seeding at the observed onset would decouple timing from growth.
5. **The elderly's attack rate.** With no elderly factor, the contact matrices alone set it, at 0.29 times the adults'. That is below PHIRST's 0.80, though PHIRST's elderly live in multigenerational households. The elderly's case counts are matched by visibility instead.
6. **Age-specific season visibility.** δ~s~ scales every age group alike, so a season visible mainly in the elderly is represented only in aggregate.

![**Figure 13. The noise budget.** Blue: each country's own week-to-week scatter around a 3-week mean, a floor no smooth curve can beat. Orange: the noise the fitted dispersion implies. The gap is misfit written off as noise.](fig13_noise_budget.png)

![**Figure 14. The seeds**, log~10~ of the infected fraction planted on 1 August. A seed ten times smaller delays a wave by about two weeks. Croatia's waves, which rise fastest, need the smallest.](fig14_seeds.png)

## 7. How far to trust it

### What the data can tell apart

| Pair | To the data | Broken by |
|:------------------------------------------|:------------------------------------------|:--------------------------------|
| R~0~ and S~0~ | exactly one number | the pin, R~0~ = 2 — judgement |
| reporting level and season visibility | exactly one number (a product) | visibility averages one |
| S~0~ and spatial spread | nearly one number: only the low tails differ | spread fixed at none — judgement |
| children's susceptibility and children's reporting | nearly one number (posterior correlation −0.31) | the cohort prior |
| season susceptibility and season visibility | different shapes (posterior correlation 0.23) | the data |
| country S~0~ and country reporting | different shapes (posterior correlation 0.10) | the data |

: **Table 6. The ties, and how each is broken.**

No direction of the likelihood is flat, and the penalised Hessian is positive definite. S~0~ is read off the rise: the observed early growth rate correlates with the fitted S~0,c,s~ at 0.46 (Spearman, 83 waves).

### Does it recover a known truth?

Data were simulated from the model onto the real design — the same countries, weeks, missing cells and populations — and refitted from scratch (Table 7, Figure 15).

| Study | The truth | What comes back |
|:-------------------------|:------------------------------|:------------------------------------------------------------|
| local, 8 replicates | the fitted optimum | country S~0~ ranking 0.97; season effects 0.97 (Pearson); visibility 0.99; 95% intervals cover 93% (median over families) |
| prior space, 4 replicates | drawn from the priors | season effects 1.00 and 0.98 (Spearman); country ranking 0.997; coverage 89–100% in every family |
| driver, 6 replicates | a covariate moves both season effects | 0.245 ± 0.029 for a true 0.25 on S~0~, and 0.359 ± 0.068 for 0.35 on visibility: unbiased |

: **Table 7. Recovery of a known truth.**

κ is not recovered, by design:

- simulated at 0.94, it comes back at 1.22 (1.10–1.36), pulled towards its prior;
- every country's S~0~ shifts about 0.06 logit down with it (the two errors correlate at −0.97), so S~0,c~ intervals cover only 77%;
- the country ranking still holds.

κ and the level of S~0~ are the cohort's; the ranking and the season effects are the data's.

![**Figure 15. Recovering a known truth** (local study, 8 replicates). Each point is one parameter in one replicate; dashed line: perfect recovery. κ (top right) comes back above its true value because its prior pulls it there.](fig15_recovery.png)

### Where it would fail

| Simulated world | Season S~0~ rank | Visibility rank | Country S~0~ rank | Noise excess |
|:----------------------------------------------------|----:|----:|----:|----:|
| control: truth the model can represent | 0.89 | 1.00 | 0.99 | 1.25× |
| true R~0~ differs between countries (sd 10%) | 0.93 | 0.93 | **0.27** | 1.27× |
| a late second wave the SIR cannot make | 0.79 | 1.00 | 0.99 | 1.27× |

: **Table 8. Misspecification tests**, one replicate each.

- **Transmissibility differences between countries:** the country ranking collapses. With R~0~ pinned, the difference has nowhere to go but S~0,c~.
- **The season results** survive both violations.
- **The noise budget** sees neither violation.

### The robustness record

Some of these checks were run on earlier versions of the model. The season effects have not moved since.

- **Evidence is judged with the wave as the unit.** Weekly residuals within a wave are autocorrelated, so raw likelihood differences, and AIC, overstate the evidence. Every comparison here uses per-wave paired differences: sign and Wilcoxon tests and bootstraps over waves, countries and seasons.
- **S~0~ or R~0~ as the sensor:** the data cannot tell a season or country effect on S~0~ from one on R~0~. The choice is a judgement.
- **The pin:** from R~0~ = 1.7 to 3.0 the fit is flat within 2 nats and the rankings are unchanged (0.99 or better); only the absolute scale moves.
- **The season asymmetry is the data's:** widening the season priors fourfold leaves both spreads unchanged. Even a deliberately hostile visibility prior leaves visibility varying nine times as much.
- **The noise gap is real:** it runs 2.1–3.7 times across 18 ways of measuring the data's own scatter.

## 8. Caveats and open questions

**Caveats for any use of the output.**

1. The country ranking of S~0~ assumes equal transmissibility and spatial structure; a 10% real difference scrambles it.
2. Absolute attack rates and reporting fractions are set by the R~0~ pin; only their patterns are data.
3. How much more susceptible children are rests on PHIRST, a South African cohort.
4. The elderly's share of infection comes from the contact matrices (0.29 times adults; PHIRST 0.80).
5. The deterministic curve does not follow within-season wander (the 2.6-fold noise gap).
6. 2025/26 is the most extreme season, on the thinnest data.
7. The bottom of the country ranking follows data quality. The three noisiest series (Poland, the Netherlands and Estonia) hold three of the four lowest S~0~ (rank correlation with dispersion 0.52). Recovery shows the estimator does not produce this pattern, so it is in the data or in misfit the simulation lacks.

**Open questions.**

- **Data.**
  - How ERVISS encodes weeks with tests but no detections (five further questions are in `documentation/to_confirm_with_surveillance.md`).
  - Whether to add Iceland, Malta and Austria.
  - Completing the 2025/26 series.
- **Model.**
  - A country-specific season deviation.
  - Seeding at the observed onset.
  - Process noise or a second strain for the noise gap.
  - A per-country R~0~, compared wave by wave, that keeps the contact matrices' own scale.
  - A European anchor for κ (Flu Watch [8]).
  - A contact-based correction for the elderly, if their attack rate becomes a target.
  - Age-specific season visibility.
- **The learning layer.** The season effects on S~0~ and on visibility are ready: recovery returns a driver's effect on both without bias. Expect little variance in susceptibility to explain (two seasons stand apart) and much in visibility. Country-level conclusions should wait for the per-country R~0~ comparison.

## References

::: {.refs}

1. Biggerstaff M, Cauchemez S, Reed C, Gambhir M, Finelli L. Estimates of the reproduction number for seasonal, pandemic, and zoonotic influenza: a systematic review of the literature. *BMC Infect Dis* 2014;14:480.
2. Cohen C, Kleynhans J, Moyes J, et al. Asymptomatic transmission and high community burden of seasonal influenza in an urban and a rural community in South Africa, 2017–18 (PHIRST): a population cohort study. *Lancet Glob Health* 2021;9:e863–74.
3. Cowling BJ, Fang VJ, Riley S, Peiris JSM, Leung GM. Estimation of the serial interval of influenza. *Epidemiology* 2009;20:344–7.
4. Gostic KM, Ambrose M, Worobey M, Lloyd-Smith JO. Potent protection against H5N1 and H7N9 influenza via childhood hemagglutinin imprinting. *Science* 2016;354:722–6.
5. Gostic KM, Bridge R, Brady S, Viboud C, Worobey M, Lloyd-Smith JO. Childhood immune imprinting to influenza A shapes birth year-specific risk during seasonal H1N1 and H3N2 epidemics. *PLoS Pathog* 2019;15:e1008109.
6. Goodwin K, Viboud C, Simonsen L. Antibody response to influenza vaccination in the elderly: a quantitative review. *Vaccine* 2006;24:1159–69.
7. Hancock K, Veguilla V, Lu X, et al. Cross-reactive antibody responses to the 2009 pandemic H1N1 influenza virus. *N Engl J Med* 2009;361:1945–52.
8. Hayward AC, Fragaszy EB, Bermingham A, et al. Comparative community burden and severity of seasonal and pandemic influenza: results of the Flu Watch cohort study. *Lancet Respir Med* 2014;2:445–54.
9. Kucharski AJ, Lessler J, Read JM, et al. Estimating the life course of influenza A(H3N2) antibody responses from cross-sectional data. *PLoS Biol* 2015;13:e1002082.
10. Monto AS, Koopman JS, Longini IM. Tecumseh study of illness. XIII. Influenza infection and disease, 1976–1981. *Am J Epidemiol* 1985;121:811–22.
11. Mossong J, Hens N, Jit M, et al. Social contacts and mixing patterns relevant to the spread of infectious diseases. *PLoS Med* 2008;5:e74.
12. Prem K, Cook AR, Jit M. Projecting social contact matrices in 152 countries using contact surveys and demographic data. *PLoS Comput Biol* 2017;13:e1005697.
13. Thompson WW, Shay DK, Weintraub E, et al. Mortality associated with influenza and respiratory syncytial virus in the United States. *JAMA* 2003;289:179–86.
14. Tsang TK, Lau LLH, Cauchemez S, Cowling BJ. Household transmission of influenza virus. *Trends Microbiol* 2016;24:123–33.

**Data.** ECDC/WHO Europe, European Respiratory Virus Surveillance Summary (ERVISS); ECDC RespiCompass. Contact matrices: Prem et al. [12]. Population and 65+ vaccination coverage: the project's data layer (`documentation/data_overview.md`).
:::

### Reproducing this report

Each step reads the previous step's output.

```
Rscript code/00_main.R                                  # model inputs from the committed data
Rscript code/07_joint_model/run_joint_model.R           # fit, diagnostics, default figures (~10 min)
Rscript code/07_joint_model/run_kappa_profile.R         # the children's modifier profile (~20 min)
Rscript code/07_joint_model/run_joint_recovery.R        # recovery and misspecification (~1.5 h)
Rscript code/07_joint_model/report/report_figures.R     # the figures of this report
sh code/07_joint_model/report/build_report.sh           # this PDF
```

The technical account is `code/07_joint_model/MODEL.md`; every decision, with its reasoning, is in `documentation/decisions.md`.
