# Data

The [project one-pager](../README_project.md#data-needs) lists the project's data needs. Each need that has
been worked on has a folder here, holding the data and a README that says where they come from, what
they cover and how far to trust them.

| Data need (one-pager) | Folder | Sources so far |
|---|---|---|
| Primary care ILI rates | [`primary_care_ili_rates/`](primary_care_ili_rates/README_primary_care_ili_rates.md) | ERVISS only, merged across its releases from 2014/15 |
| Primary care test positivity for flu | [`primary_care_flu_positivity/`](primary_care_flu_positivity/README_primary_care_flu_positivity.md) | ERVISS only, merged across its releases from 2014/15 |
| SARI rates | [`sari_rates/`](sari_rates/README_sari_rates.md) | ERVISS only (merged from 2020-W40, its first SARI week, in the hospitalisation folder) |
| SARI test positivity for flu | [`sari_flu_positivity/`](sari_flu_positivity/README_sari_flu_positivity.md) | ERVISS only (merged from 2020-W40, its first SARI week, in the hospitalisation folder) |
| Diverse influenza hospitalisation indicators | [`flu_hospitalisation_indicators/`](flu_hospitalisation_indicators/README_flu_hospitalisation_indicators.md) | ERVISS, merged across releases, plus seven national sources from a wider search |
| Dominant influenza subtype | not yet | an earlier analysis's table is in `external/` (below) |
| Influenza vaccination coverage | not yet | RespiCompass and an external 65+ table (below) |

## Other files in `data/`

These are inputs of the analysis pipeline that do not yet have a data-need folder.
[`documentation/data_overview.md`](../documentation/data_overview.md) describes how the pipeline loads
them.

| File | Contents |
|---|---|
| `slim_flu_iliplus.csv` | the combined influenza ILI+ panel (ILI rate × positivity), built by [`build_slim_panel.R`](../code/04_modelling/build_slim_panel.R) with the rules in [`stitch_iliplus.R`](../code/01_main_supporting/stitch_iliplus.R) |
| `data_respicompass_iliplus.csv` | RespiCompass's influenza ILI+, the panel's source up to 2023/24 |
| `erviss_flu_type_subtype.csv` | ERVISS influenza type and subtype breakdown; read by [`dominant_subtype.R`](../code/05_analysis/dominant_subtype.R) |
| `erviss_sequencing.csv`, `erviss_variants.csv` | ERVISS sequencing volumes and variant shares |
| `external/` | drivers sourced outside ERVISS: dominant subtype by season, vaccination coverage 65+ after COVID-19, vaccine effectiveness, winter climate; provenance in [`external_drivers.md`](../documentation/external_drivers.md) |
| `vax_flu_history.csv`, `vax_flu_history_all.csv`, `vax_flu_scenarios.csv` | RespiCompass vaccination coverage: history (65+, and all target groups) and scenarios |
| `MUestimates_all_locations_{1,2}.xlsx` | synthetic contact matrices (Prem et al.) |
| `population_pyramid.fst` | population by age (ECDC snapshot) |
| `vax_flu_data.csv`, `survey_deidentified.xlsx`, `from Nick end Nov test respicast/` | not read by the code |
