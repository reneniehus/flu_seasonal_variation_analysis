# Primary care ILI rates

Data need from the [project one-pager](../../README.md#data-needs): *primary care ILI rates
[season, country, week, age]*.

**Source.** ERVISS (ECDC), file
[`ILIARIRates.csv`](https://github.com/EU-ECDC/Respiratory_viruses_weekly_data/tree/main/data). It
holds weekly consultation rates for influenza-like illness (ILI) and acute respiratory infection (ARI)
from sentinel primary-care networks, by age group.

**Search scope.** Only ERVISS was used. No wider search has been done for this data need.

## The file here

[`erviss_iliari.csv`](erviss_iliari.csv) is the analysis pipeline's snapshot of ERVISS
`ILIARIRates.csv`. It is unchanged from ERVISS, committed on 2026-06-12, and covers 2021-W25 to
2026-W20.

- **Columns:** `survtype, countryname, yearweek, indicator, age, value`.
- **Indicators:** `ILIconsultationrate` and `ARIconsultationrate`.
- **Ages:** `0-4, 5-14, 15-64, 65+, total`.
- **Loading:** [`load_data.R`](../../code/01_main_supporting/load_data.R) loads it as
  `data$epi$erviss_ili_ari`, and `00_main.R` reads it from this file (`new_from_online = F`). Refreshing
  it from ERVISS (`new_from_online = T`) would replace it with ERVISS's current rolling window, which
  starts at 2022-W25 and so would drop 2021/22.

**Units.** ERVISS reports rates per 100 000 population, with exceptions:
- Cyprus, Luxembourg and Malta: per 100 consultations;
- Finland: per 100 000 consultations.

The ILI+ panel multiplies Cyprus, Luxembourg and Malta by 1000. Whether that is right is
[open question 6](../../documentation/to_confirm_with_surveillance.md).

## Coverage

ILI rate, all ages: weeks with data per season.

| Country | 2021/22 | 2022/23 | 2023/24 | 2024/25 | 2025/26 |
|---|---:|---:|---:|---:|---:|
| Austria | 27 | 26 | 27 | 27 | 27 |
| Belgium | 52 | 52 | 52 | 52 | 33 |
| Croatia | 37 | 36 | 51 | 50 | 33 |
| Czechia | 52 | 52 | 52 | 52 | 33 |
| Denmark | 52 | 52 | 52 | 52 | 33 |
| Estonia | 52 | 51 | 46 | 49 | 30 |
| Finland | 52 | 51 | 49 | 51 | 28 |
| France | 41 | 39 | 38 | 37 | 26 |
| Greece | 52 | 52 | 52 | 52 | 33 |
| Hungary | 33 | 33 | 32 | 32 | 32 |
| Iceland | 33 | 51 | 47 | 51 | 32 |
| Ireland | 52 | 51 | 52 | 52 | 33 |
| Italy | 28 | 28 | 27 | 28 | |
| Latvia | 36 | 32 | 32 | 33 | 32 |
| Lithuania | 52 | 51 | 52 | 52 | 33 |
| Luxembourg | 50 | 52 | 52 | 52 | 33 |
| Malta | 49 | 51 | 33 | 52 | 33 |
| Netherlands | 52 | 52 | 52 | 52 | |
| Norway | 43 | 52 | | 52 | 33 |
| Poland | 50 | 52 | 52 | 52 | 33 |
| Portugal | 52 | 21 | | | |
| Romania | 51 | 48 | 51 | 50 | 32 |
| Slovakia | 52 | 52 | 52 | 52 | 33 |
| Slovenia | 52 | 52 | 52 | 52 | 33 |
| Spain | | | 52 | 52 | 33 |

**Who reports what.**
- **25 countries** report ILI rates. All except Greece give them by age. 21 also report ARI rates.
- **No ILI rates:** Bulgaria, Germany and Sweden do not report ILI rates to ERVISS. Italy and the
  Netherlands stop after 2024/25.
- **Winter only:** Austria, Hungary, Italy and Latvia report 26–36 weeks a season.
- **Incomplete 2025/26:** it ends at 2026-W20 in this snapshot
  ([open question 5](../../documentation/to_confirm_with_surveillance.md)).

## Earlier seasons

This file starts at 2021-W25. For 2014/15–2023/24 the project's ILI+ panel uses RespiCompass's ILI+,
[`data_respicompass_iliplus.csv`](../data_respicompass_iliplus.csv), in which ILI × positivity is
already combined. The stitching rules are in
[`stitch_iliplus.R`](../../code/01_main_supporting/stitch_iliplus.R), and the resulting panel is
[`slim_flu_iliplus.csv`](../slim_flu_iliplus.csv).

**A lead within ERVISS, not yet used.** ERVISS's own dated snapshots from 2024-01-26 to 2024-10-18
carry ILI rates back to 2014-W40 for 28 countries, 23 of them by age group. Merging them, as was done
for the [hospital indicators](../flu_hospitalisation_indicators/README.md#1-erviss), would give
ERVISS's own pre-COVID-19 seasons. That would check RespiCompass's ILI+ and add age groups.

## Related

- [Primary care test positivity](../primary_care_flu_positivity/README.md): the other half of ILI+.
- [Open questions about these data](../../documentation/to_confirm_with_surveillance.md).
