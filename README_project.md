---
editor_options: 
  markdown: 
    wrap: 72
---

*Project on a page*

# Analysing seasonal variation in influenza burden across the EU/EEA

**Date:** October 2026

## Project team

-   **Rene** — Project lead; runs the analysis, writes outputs
-   **Edoardo** — Line manager, supervisor and internal end-user
-   **Nick, Luisa, Leah** — Surveillance data interpretation and
    context. **Input needed**
-   **VPI team** — VE estimates (VPI team)

## Justification (problem/opportunity, alignment with team and organisation goals)

-   Before each winter, EU/EEA authorities plan vaccination and
    healthcare capacities in response to seasonal virus burden. This
    planning has to be done in the absence of reliable burden
    predictions.
-   This is because forecasts, such as those provided by ECDC's
    RespiCast, are only short term, and long-term predictions of viral
    burden has proved highly challenging.
-   As a step toward better understanding viral circulation and seasonal
    burden patterns, we focus here on influenza and carefully analyse
    historical variation of relevant indicators across seasons and
    across the EU/EEA.
-   By examining the association between such variation and plausible
    drivers of viral burden, we generate hypothesis about the
    predictable nature of seasonal influenza burden.

## Data needs

-   [Primary care ILI
    rates](data/primary_care_ili_rates/README_primary_care_ili_rates.md)
    `[season,country,week,age]`
-   [Primary care test positivity for
    flu](data/primary_care_flu_positivity/README_primary_care_flu_positivity.md)
    `[season,country,week]`
-   [SARI rates](data/sari_rates/README_sari_rates.md)
    `[season,country,week,age]`
-   [SARI test positivity for
    flu](data/sari_flu_positivity/README_sari_flu_positivity.md)
    `[season,country,week]`
-   [Diverse influenza hospitalisation
    indicators](data/flu_hospitalisation_indicators/README_flu_hospitalisation_indicators.md)
    `[season,country]`
-   Dominant influenza subtype `[season,country]`
-   Influenza vaccination coverage `[season,country]`

## Objectives (SMART: specific, measurable, achievable, relevant, time-bound) and Deliverables

1.  Determine relevant influenza indicators: we establish for each
    EU/EEA country (where possible) a surveillance indicator that
    correlates with severe disease burden caused by influenza. This is
    done through combining evidence from temporal correlation between
    age-disaggregated indicators (including all-cause mortality,
    secondary care, and primary care indicators) and discussions with
    surveillance and disease experts. The result is a 1-page report by
    mid-January 2027 on\
    surveillance indicators that correlate with severe disease burden
    caused by influenza.
2.  Describing patterns and trends of influenza indicators: we describe
    dynamics of influenza indicators across the EU/EEA using
    season-specific attributes (e.g., onset timing, peak timing,
    steepness). We then examine how those attributes relate to one
    another across seasons and across the EU/EEA, we examine how
    specific seasons were different across the EU/EEA, and how specific
    countries were different to others. We perform this analysis without
    imposing any mechanistic representation of viral spread (e.g., no
    mechanistic model fitting). The result is a 2-page descriptive
    report of influenza indicator dynamics for all available
    country-seasons.
3.  Correlation analysis and building causal hypotheses: we correlate
    the observed patterns of influenza burden indicators with plausible
    biological drivers and descriptors or surveillance systems. Through
    discussing the results with disease experts and re-adjusted analyses
    of correlations, we generate hypothesis about the mechanistic
    drivers underlying influenza burden indicator. This yields
    discussions of how these mechanisms could be further investigated
    and how they could aid predictability of future influenza burden.
    The result is a 2-4 page technical brief summarising methods,
    results, discussion and future steps.

## Activity scope

**In scope**

-   Review of existing literature
-   Analysis of ERVISS data outputs and relevant TESSY data
-   Scraping of EU/EEA national sources on influenza hospitalisation
    data (including validation through discussion with surveillance data
    expert)
-   Discussions with surveillance data expert colleagues
-   Phenomenological descriptors of seasonal dynamics, with applied
    smoothing where appropriate
-   Building causal hypothesis from correlations
-   Discussion of caveats and potential confounders

**Out of scope**

-   Systematic literature review
-   Cleaning of data streams where ERVISS alternatives exist
-   Engaging with Network representatives
-   Mechanistic model building or fitting
-   New threshold analyses or trend analyses
-   Scenario projections
-   Analysis of COVID-19, RSV or other disease indicators

## Dependencies / constraints

-   **Timely input from named colleagues**\
    Responsible: Project lead; named colleagues
-   **Lack of influenza-confirmed hospitalisation data in the EU/EEA;
    cross-country surveillance differences; confounds pre/post-COVID**\
    Responsible: Project lead
-   **Limited capacity of the project lead (currently \<0.5 FTE for
    modelling)**\
    Responsible: Project lead; line manager

*Classified as ECDC NORMAL*

------------------------------------------------------------------------

**In this repository.**

-   **Data:** each data need above links to its folder under
    [`data/`](data/README_data.md).
-   **Code and methods, and how to run them:**
    [`code/README_code.md`](code/README_code.md).
-   **Longer background and references:**
    [`documentation/PROJECT_SCOPE.md`](documentation/PROJECT_SCOPE.md).
-   **Why things were done as they were:**
    [`documentation/decisions.md`](documentation/decisions.md).
