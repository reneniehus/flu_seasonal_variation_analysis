# Drivers of inter-seasonal variability in seasonal influenza across the EU/EEA: a retrospective analysis

> The project is summarised on one page in the [repository README](../README.md). This document
> gives the longer background, the research questions and the references.

## Background and rationale

Seasonal influenza recurs every winter across the EU/EEA, yet the shape of each epidemic varies markedly — both between countries and from one year to the next — in the timing of onset and peak, peak intensity, duration, and cumulative burden. Each year, well before a season begins, public health authorities across the EU/EEA must commit to preparations: communicating the importance of vaccination, procuring and allocating vaccines, and planning hospital and ICU capacity. These commitments are made while the course of the coming season is still largely unknown. Preparing under this uncertainty is costly: it risks misallocating scarce healthcare resources, or leaving communities across the EU/EEA under-protected when a severe season arrives — the citizens whose health ECDC is mandated to protect.

Short-term forecasting is an active field of research aimed at improving the capability of model forecasts to capture observed influenza indicators a few weeks ahead. ECDC's RespiCast — the European Respiratory Diseases Forecasting Hub — and the US CDC's FluSight challenge both aggregate predictions from multiple modelling teams into ensembles, typically at horizons of one to four weeks [1, 2]. Such ensembles consistently outperform individual models, but forecast uncertainty is large and increases rapidly toward the four-week horizon [1, 2]. These methods strengthen within-season situational awareness, but they do not provide the months-ahead foresight that procurement and capacity decisions require.

Scenario and medium-term modelling — for example through RespiCompass, ECDC's respiratory-disease scenario modelling hub — explores plausible trajectories under defined assumptions. Multi-model exercises, however, reveal substantial disagreement between models arising from differing structural assumptions about transmission and immunity [3, 4]. More fundamentally, the parameters that most shape a season — which (sub)type dominates, how closely the vaccine antigenically matches circulating strains, and real-world vaccine effectiveness — are not knowable far in advance. This is an intrinsic limit on long-range prediction, not merely a modelling shortcoming.

A useful first step is therefore to look backwards. Rather than predict the next season, this project asks what has explained the variability already observed across EU/EEA seasons. A substantial literature links inter-seasonal differences to a recurring set of factors: the dominant type/subtype, with A(H3N2)-dominated seasons associated with greater morbidity and mortality, particularly in older adults [5]; antigenic novelty and vaccine match, with H3N2 seasons tending to show lower vaccine effectiveness through more frequent mismatch [6, 7]; vaccine coverage [4]; prior-season immunity and population susceptibility — including the post-pandemic disruption to epidemic timing and intensity that followed the relaxation of non-pharmaceutical interventions [8]; demography and contact patterns; and climatic conditions such as temperature and humidity [9]. Recent work demonstrates that these drivers can be systematically related to epidemic burden, timing, and the age distribution of cases within a single setting [7]. An equivalent, EU/EEA-wide characterisation is currently lacking — and is a prerequisite for deciding which factors deserve a place in any future predictive effort.

## Aim and research questions

The aim is to characterise the inter-seasonal and inter-country variability of seasonal influenza epidemics in the EU/EEA, and to identify candidate factors that explain it.

- How have past EU/EEA influenza seasons varied — across countries and years — in onset and peak timing, peak intensity, duration, and cumulative burden?
- Which candidate drivers (dominant (sub)type, antigenic match and vaccine effectiveness, vaccine coverage, prior-season activity, demography, climate) are associated with these seasonal features?
- How could the associations identified be used to better understand, anticipate, and prepare for future influenza seasons in the EU/EEA?

The analysis is exploratory and hypothesis-generating: it seeks robust, plausible associations rather than causal estimates.

## Approach

The strategy is to turn each country–season influenza epidemic into a small set of interpretable
**summary features**, and then explore how those features relate to candidate drivers.

1. **Characterise the epidemics (RQ1).** From ECDC ERVISS sentinel surveillance we build weekly
   **ILI+** curves (influenza-like-illness consultation rate × influenza positivity) per country and
   season, and extract features describing each wave: onset and peak timing, peak intensity,
   duration / cumulative burden, and growth rate. Features are produced by a small set of
   **complementary methods** that share one output schema, so each feature can be cross-checked
   rather than trusted from a single model:
   - a **susceptible-reconstruction SIR**, fitted both deterministically and as an Extended Kalman
     Filter, summarising a season by an initial **susceptibility** `S0` (read off the wave's rise
     rate) and a reporting fraction;
   - a **descriptive** method that smooths the curve and reads timing, intensity and steepness
     directly as observed-shape features (no SIR fit, no susceptibility mapping).
   The SIR here is a *summarising* device — a susceptibility reconstruction — not a causal
   transmission model (see Out of scope): transmission parameters are fixed from the literature and
   only the *relative* susceptibility across seasons is interpreted.

2. **Assemble candidate drivers (RQ2).** For each country–season we collate plausible drivers —
   dominant (sub)type, antigenic match / vaccine effectiveness, vaccine coverage, prior-season
   activity, demography, and climate — from ERVISS / RespiCompass and external sources.

3. **Explore associations (RQ2–3).** We relate the epidemic features to the drivers across countries
   and seasons, visually and with simple, transparent statistics (rank correlations / regressions),
   reporting robust, plausible associations and interpreting them for preparedness — i.e. which
   factors most deserve attention when anticipating a season.

*Current state:* the surveillance data layer, the feature-extraction methods, and a first pass of the
within-country driver analysis are implemented (`code/05_analysis/`; results in
`documentation/findings_descriptors.md`, reasoning in `documentation/analysis_strategy.md`). Further
drivers (climate, subtype novelty, vaccine effectiveness) are the next stage.

## In scope

- Retrospective characterisation of past EU/EEA influenza seasons by per-country, per-season summary
  features (onset / peak timing, intensity, duration, burden, susceptibility), extracted with
  complementary methods.
- Exploratory association of those features with plausible drivers (dominant (sub)type, vaccine
  coverage, antigenic match and vaccine effectiveness, prior-season immunity, demography, climate).
- Hypothesis generation on which factors drive influenza epidemiology and should therefore inform its
  prediction.
- A structured literature review of related work.

## Out of scope

- Influenza seasons outside the EU/EEA, and Southern Hemisphere dynamics (the constant seeding
  assumption stands in for off-season importation without modelling it).
- Operational forecasts and long-term projections.
- Causal or mechanistic transmission modelling — the analysis is associational and
  hypothesis-generating by design; the SIR is used only to summarise each curve into an interpretable
  susceptibility.

## Repository

How to run, the data layer and the method framework are documented in `code/README.md` and
`documentation/` (quickstart, data overview, the model maths). The data needs, each with its sources
and coverage, are under `data/` (index: `data/README.md`). Key design choices and their rationale are
recorded in `documentation/decisions.md`.

## Key references and related work

1. ECDC. *RespiCast — the European Respiratory Diseases Forecasting Hub* (ensemble forecasts of ILI/ARI incidence, 1–4 weeks ahead; 2023/24–2024/25 seasons). https://respicast.ecdc.europa.eu/
2. Reich NG, et al. A collaborative multiyear, multimodel assessment of seasonal influenza forecasting in the United States (FluSight). *PNAS*. 2019;116(8):3146–54. https://www.pnas.org/doi/10.1073/pnas.1812594116 — current hub: CDC FluSight Forecast Hub, https://github.com/cdcepi/FluSight-forecast-hub
3. ECDC. *RespiCompass — European respiratory diseases scenario modelling hub.* https://respicompass.ecdc.europa.eu/
4. ECDC. *Multi-model analysis to quantify the impact of vaccination on COVID-19 and influenza hospitalisation burden among older adults in the EU/EEA, 2024/25.* Stockholm: ECDC; 2025. doi:10.2900/3212996
5. *Subtype H3N2 Influenza A Viruses: An Unmet Challenge in the Western Pacific.* Vaccines. 2022;10(1):112. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8778411/
6. *Multi-strain modeling of influenza vaccine effectiveness in older adults and its dependence on antigenic distance.* Scientific Reports. 2024. https://www.nature.com/articles/s41598-024-72716-1
7. Perofsky AC, et al. *Antigenic drift and subtype interference shape A(H3N2) epidemic dynamics in the United States.* eLife. 2024. https://elifesciences.org/reviewed-preprints/91849
8. *Global analysis of influenza epidemic characteristics in the first two seasons after lifting non-pharmaceutical interventions for COVID-19.* International Journal of Infectious Diseases. 2024. https://www.ijidonline.com/article/S1201-9712(24)00447-8/fulltext
9. Tamerius J, et al. *Environmental Predictors of Seasonal Influenza Epidemics across Temperate and Tropical Climates.* PLoS Pathogens. 2013;9(3):e1003194. https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1003194
