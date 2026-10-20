# Data overview

What `load_data()` assembles into the `data` list, and what `gen_model_input()` turns it into.
All streams are weekly, EU/EEA countries, keyed by ISO2 country code (`country_short`) and an
ISO-week Wednesday `date`. Seasons run Aug 1 -> Jul 31 (configurable in the settings).

## `data` — raw, per-source streams

### `data$epi` — ECDC ERVISS + RespiCompass surveillance
Loaded via the registry in `code/01_main_supporting/load_data.R` (add a row -> new stream).
Source: <https://github.com/EU-ECDC/Respiratory_viruses_weekly_data>. Snapshots that serve a data need
of the project one-pager sit in that need's folder under `data/`, with a README on source, coverage
and caveats (index: `data/README_data.md`).

| key | file | schema | contents |
|---|---|---|---|
| `erviss_ili_ari` | ILIARIRates.csv | rates | ILI & ARI primary-care consultation rates, by age |
| `erviss_sari_rates` | SARIRates.csv | rates | SARI hospitalisation rates, by age |
| `erviss_typing_sentinel` | sentinelTestsDetectionsPositivity.csv | detailed | sentinel tests/detections/positivity per pathogen |
| `erviss_typing_nonsentinel` | nonSentinelTestsDetections.csv | detailed | non-sentinel tests/detections per pathogen |
| `erviss_typing_sari` | SARITestsDetectionsPositivity.csv | detailed | SARI tests/detections/positivity per pathogen |
| `erviss_flu_type_subtype` | activityFluTypeSubtype.csv | detailed | influenza type/subtype breakdown |
| `erviss_severity_nonsentinel` | nonSentinelSeverity.csv | detailed | severity indicators (non-sentinel) |
| `erviss_sequencing` | sequencingVolumeDetectablePrevalence.csv | detailed | sequencing volume / detectable prevalence |
| `erviss_variants` | variants.csv | detailed | variant shares |
| `respicompass_iliplus` | RespiCompass ili_plus.csv | rates | externally-provided influenza ILI+ |

- **rates** schema = slim columns: `country_short, date, target, agegroup, value`.
- **detailed** schema = all original columns kept (`pathogen, pathogentype, pathogensubtype,
  indicator, age, value, ...`) plus added `date` and `country_short`.
- Pathogens covered by the typing streams: **Influenza, SARS-CoV-2, RSV**.

### Other `data$` slots
| slot | contents | source |
|---|---|---|
| `data$vax` | `data_vax` (forward coverage scenarios), `data_vax_history` (65+), `data_vax_history_all` (all target groups) | RespiCompass |
| `data$contact` | per-country synthetic contact matrices (list) | Prem et al. (shipped in `data/`) |
| `data$helpers_respicompass` | `iso2_code`, `iso_weeks` lookup tables | RespiCompass |
| `data$demography_respicast` | `population_pyramid`, `population_pyramid_fine` | RespiCompass |
| `data$demography_ECDC` | `population_pyramid` | committed snapshot `data/population_pyramid.fst` |

## `models_in` — tidy, model/plot-ready tables (`gen_model_input()`)

- **`data_timeseries_long`** — the single source of truth. One row per
  `country_short × season × date × source × stream × pathogen × indicator × agegroup`, with
  `value`, `unit`, `observed`, and season-time columns (`season_week`, `iso_week`, ...).
  Indicators: `ILIconsultationrate`, `ARIconsultationrate`, `tests`, `detections`,
  `positivity`, `ili_plus` (per pathogen), `vaccine_coverage`.
  Streams include `ili_ari`, `typing_sentinel/nonsentinel`, `ili_plus_sentinel/nonsentinel/
  respicompass`, `vaccination_history_65plus`, `vaccination_scenario`.
- **`data_timeseries_wide`** — weekly streams pivoted to one column per indicator series,
  keyed by `country_short × season × date × agegroup`.
- **`data_season_summary`** — per `country × season × series` quality/summary stats:
  `n_weeks_observed`, `weeks_in_season`, `weeks_in_span`, `completeness`,
  `completeness_active`, `sum/mean/max_value`, `peak_date`, `first/last_date`
  (age-specific and pooled-over-age rows).
- **`contacts`** — transformed 4-age-group contact matrices per country (+ EU average).

### Derived "+" indicators
`ili_plus` = ILI consultation rate × pathogen positivity, built for Influenza / SARS-CoV-2 /
RSV from both sentinel and non-sentinel typing. Positivity is `detections / tests` at the
pathogen-total level (`pathogentype == pathogen`); it is age-total, broadcast across the
age-specific ILI rate.

### Two completeness measures (see captions in the eyeballing report)
- `completeness` — reported weeks / ISO weeks in the **full Aug–Jul season** (penalises the
  off-season; winter-only streams cap well below 100%).
- `completeness_active` — reported weeks / weeks between the **first and last reported week**
  (off-season ignored; isolates genuine mid-season gaps).
