# SARI rates

Data need from the [project one-pager](../../README.md#data-needs): *SARI rates [season, country, week,
age]*.

**Source.** ERVISS (ECDC), file `SARIRates.csv`. It holds weekly rates of severe acute respiratory
infection (SARI) among patients admitted to sentinel hospitals, by age group. The rate counts SARI of
any cause, so on its own it is not influenza-specific. Multiplied by the influenza positivity of SARI
samples it becomes SARI+; see [SARI test positivity](../sari_flu_positivity/README.md) and the
[hospital indicators](../flu_hospitalisation_indicators/README.md).

**Search scope.** Only ERVISS was used. No wider search has been done for this data need.

## The file here

[`erviss_sari_rates.csv`](erviss_sari_rates.csv) is the analysis pipeline's snapshot of ERVISS
`SARIRates.csv`. It is unchanged from ERVISS, committed on 2026-06-12, and covers 2021-W25 to
2026-W20.

- **Columns:** `survtype, countryname, yearweek, indicator, age, value`.
- **Indicator:** `SARIrate`.
- **Ages:** `0-4, 5-14, 15-64, 65+, total`.
- **Loading:** [`load_data.R`](../../code/01_main_supporting/load_data.R) loads it as
  `data$epi$erviss_sari_rates`.

**Units.** The rate is per 100 000 population. Countries estimate the population their sentinel
hospitals serve in different ways, so levels compare within a country, not between countries.

## Coverage

SARI rate, all ages: weeks with data per season.

| Country | 2021/22 | 2022/23 | 2023/24 | 2024/25 | 2025/26 |
|---|---:|---:|---:|---:|---:|
| Austria | | 21 | 52 | 52 | 33 |
| Belgium | 52 | 52 | 52 | 52 | 32 |
| Croatia | 39 | 52 | 39 | | |
| Cyprus | | | 1 | 52 | 33 |
| Estonia | 52 | 52 | 52 | 52 | 32 |
| Germany | 52 | 52 | 52 | 52 | 32 |
| Greece | | 44 | 52 | 52 | 33 |
| Iceland | | | | 38 | 17 |
| Ireland | 51 | 52 | 52 | 52 | 33 |
| Latvia | | | | 23 | 32 |
| Lithuania | | 27 | 47 | | |
| Luxembourg | | 40 | 52 | 52 | 33 |
| Malta | 52 | 52 | 52 | 52 | 33 |
| Romania | 26 | 27 | 27 | 27 | 27 |
| Slovakia | | 51 | 52 | 52 | 32 |
| Spain | 52 | 52 | 52 | 52 | 33 |

**Who reports what.**
- **16 countries** report a SARI rate, 14 of them by age (all but Cyprus and Greece).
- **Austria** reports the rate but no SARI testing, so it has no SARI+.
- **Germany's** rate is RKI's ICOSARI series of all SARI cases.
- **Romania** reports winter weeks only.

## A longer series

**Earliest SARI data.** ERVISS's SARI data start at 2020-W40. The
[hospital indicators](../flu_hospitalisation_indicators/README.md) folder merges every SARI snapshot
since 2024-11 with the current file. Its
[`erviss_flu_hospital_weekly.csv`](../flu_hospitalisation_indicators/erviss_flu_hospital_weekly.csv)
(indicator `sari_rate`) therefore runs from 2020-W40 to 2026-W40, longer than the snapshot here.

**National sources reaching further back.**
- Germany (RKI ICOSARI): the all-cause SARI rate from 2014/15.
- Spain (SiVIRA): from 2020/21.

Both were found by the hospital-indicator search; neither is merged here.

## Related

- [SARI test positivity](../sari_flu_positivity/README.md): makes this rate influenza-specific
  (SARI+).
- [Influenza hospitalisation indicators](../flu_hospitalisation_indicators/README.md): SARI+ for each
  country, with the national series.
