*Project on a page*

# Analysing seasonal variation in influenza burden across the EU/EEA

**Date:** October 2026

## Project team

- **Rene** — Project lead; runs the analysis, writes outputs
- **Edoardo** — Line manager, supervisor and internal end-user
- **Nick, Luisa, Leah** — Surveillance data interpretation and context. **Input needed**
- **VPI team** — VE estimates (VPI team)

## Justification (problem/opportunity, alignment with team and organisation goals)

- Before each winter, EU/EEA authorities plan vaccination and healthcare capacities in response to
  seasonal virus burden in the absence of reliable predictions.
- Forecasts provided by ECDC's RespiCast are only short term, and long-term prediction of viral burden
  has proved highly challenging
- As a step towards understanding viral circulation and burden patterns, we focus here on influenza
  and carefully analyse historical variation of relevant indicators across seasons and across the
  EU/EEA
- By examining the association between such variation and plausible drivers of viral burden, we
  generate hypothesis about the predictable nature of seasonal influenza burden.

## Data needs

- [Primary care ILI rates](data/primary_care_ili_rates/README_primary_care_ili_rates.md) \[season,country,week,age\]
- [Primary care test positivity for flu](data/primary_care_flu_positivity/README_primary_care_flu_positivity.md) \[season,country,week\]
- [SARI rates](data/sari_rates/README_sari_rates.md) \[season,country,week,age\]
- [SARI test positivity for flu](data/sari_flu_positivity/README_sari_flu_positivity.md) \[season,country,week\]
- [Diverse influenza hospitalisation indicators](data/flu_hospitalisation_indicators/README_flu_hospitalisation_indicators.md) \[season,country\]
- Dominant influenza subtype \[season,country\]
- Influenza vaccination coverage \[season,country\]

## Objectives (SMART: specific, measurable, achievable, relevant, time-bound) and Deliverables

1. Determine surveillance indicators that correlate with severe disease burden caused by influenza: by
   mid-December 2026 produce a 1-page report that establishes a selection of influenza indicators
   based on discussion with surveillance and disease experts and evidence on temporal correlation
   between influenza-related indicators.
2. By end of February 2027, produce a descriptive 2-page report of relevant influenza indicator(s)
   dynamics for all country-seasons.
3. By end of March 2027, examine whether plausible seasonal drivers of influenza burden relate to
   between season variation. Through discussing these results with disease experts, this may generate
   hypothesis about the predictability of future influenza season.
4. The resulting hypotheses and supporting evidence are presented in a 2–4-page technical brief.

## Scope

**In scope**

- Influenza ILI+ (ILI rate × positivity), 2014/15–2018/19 and 2023/24–2025/26; panel refreshed to
  complete 2025/26
- Burden compared within country (reporting scales differ by orders of magnitude); timing and growth
  also across countries
- Biphasic-season flags; robustness to smoothing, onset threshold and RespiCompass-to-ERVISS source
  change (2023/24 the only overlap season)
- Indicator screen with existing within-country models

**Out of scope**

- Forecasts or scenario projections
- Causal or mechanistic modelling
- Analysis of COVID-19, RSV or other disease indicators

## Dependencies / constraints

- **Timely input from named colleagues**\
  Responsible: Project lead; named colleagues
- **Lack of influenza-confirmed hospitalisation data in the EU/EEA; cross-country surveillance
  differences; confounds pre/post-COVID**\
  Responsible: Project lead
- **Limited capacity of the project lead (currently <0.5 FTE for modelling)**\
  Responsible: Project lead; line manager

*Classified as ECDC NORMAL*

---

**In this repository.**
- **Data:** each data need above links to its folder under [`data/`](data/README_data.md).
- **Code and methods, and how to run them:** [`code/README_code.md`](code/README_code.md).
- **Longer background and references:**
  [`documentation/PROJECT_SCOPE.md`](documentation/PROJECT_SCOPE.md).
- **Why things were done as they were:** [`documentation/decisions.md`](documentation/decisions.md).
