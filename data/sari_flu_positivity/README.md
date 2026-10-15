# SARI test positivity for influenza

Data need from the [project one-pager](../../README.md#data-needs): *SARI test positivity for flu
[season, country, week]*.

**Source.** ERVISS (ECDC), file `SARITestsDetectionsPositivity.csv`. It holds weekly influenza tests,
detections and positivity among patients with severe acute respiratory infection (SARI) in sentinel
hospitals, by type and subtype, all ages.

**Search scope.** Only ERVISS was used. No wider search has been done for this data need.

## The file here

[`erviss_detections_sari.csv`](erviss_detections_sari.csv) is the analysis pipeline's snapshot of
ERVISS `SARITestsDetectionsPositivity.csv`. It is unchanged from ERVISS, committed on 2026-06-12, and
covers 2021-W25 to 2026-W20.

- **Columns:** `survtype, countryname, yearweek, pathogen, pathogentype, pathogensubtype, indicator,
  age, value`.
- **Indicators:** `tests`, `detections`, `positivity`.
- **Pathogens:** influenza, SARS-CoV-2 and RSV.
- **Loading:** [`load_data.R`](../../code/01_main_supporting/load_data.R) loads it as
  `data$epi$erviss_typing_sari`.

## Coverage

Influenza tests among SARI patients: weeks reported per season, and the median number of tests in a
winter week (weeks 40–20).

| Country | 2021/22 | 2022/23 | 2023/24 | 2024/25 | 2025/26 | Median tests a winter week |
|---|---:|---:|---:|---:|---:|---:|
| Belgium | 52 | 52 | 52 | 52 | 32 | 138 |
| Croatia | 39 | 52 | 39 | | | 30 |
| Cyprus | | | 1 | 52 | 33 | 57 |
| Estonia | 52 | 52 | 52 | 52 | 33 | 122 |
| Germany | | 52 | 52 | 52 | 32 | (690, not real tests: see below) |
| Greece | | 44 | 52 | 52 | 33 | 53 |
| Iceland | | | | 39 | 17 | 7 |
| Ireland | 39 | 52 | 52 | 52 | 33 | 22 |
| Latvia | | | | 23 | 32 | 481 |
| Lithuania | | 27 | 47 | | | 36 |
| Luxembourg | | | | 39 | 27 | 24 |
| Malta | 52 | 52 | 52 | 52 | 33 | 38 |
| Romania | 27 | 29 | 28 | 27 | 29 | 21 |
| Slovakia | | 23 | 24 | 23 | 15 | 3 |
| Spain | | 52 | 52 | 52 | 33 | 375 |

**Caveats by country.**
- **15 countries** report SARI testing, plus an EU/EEA aggregate.
- **Germany:** the "positivity" is not a laboratory result. It is the share of SARI cases coded ICD-10
  J10 (influenza, virus identified), and its "tests" are a pseudo-count. This was established in the
  hospital-indicator verification.
- **Austria** reports SARI rates but no testing.
- **Few tests:** Slovakia (about 3 a week) and Iceland (about 7) test too few patients for a stable
  weekly positivity.
- **Ages:** positivity is published for all ages only. Applying it to age-specific SARI rates assumes
  it is equal across ages (owner decision, 2026-10).

## Where it is used

SARI+ = SARI rate × influenza positivity / 100. It is built from merged ERVISS releases (2020-W40 to
2026-W40, longer than this snapshot) in the
[hospital indicators](../flu_hospitalisation_indicators/README.md) folder. There
[`erviss_flu_hospital_weekly.csv`](../flu_hospitalisation_indicators/erviss_flu_hospital_weekly.csv)
has indicators `sari_flu_tests`, `sari_flu_detections` (also by subtype), `sari_flu_positivity` and
`sari_plus`.

## Related

- [SARI rates](../sari_rates/README.md): the other half of SARI+.
- [Primary care test positivity](../primary_care_flu_positivity/README.md): the same measure among
  primary-care patients.
