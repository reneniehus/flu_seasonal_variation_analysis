# Primary care test positivity for influenza

Data need from the [project one-pager](../../README.md#data-needs): *primary care test positivity for
flu [season, country, week]*.

**Source.** ERVISS (ECDC), file `sentinelTestsDetectionsPositivity.csv`. It holds weekly influenza
tests, detections and positivity for swabs taken from patients seen in sentinel primary-care networks,
by type and subtype.

**Search scope.** Only ERVISS was used. No wider search has been done for this data need.

## The files here

Both files are the analysis pipeline's snapshots: unchanged from ERVISS, committed on 2026-06-12, and
covering 2021-W25 to 2026-W20. Both list influenza, SARS-CoV-2 and RSV.

- **Columns:** `survtype, countryname, yearweek, pathogen, pathogentype, pathogensubtype, indicator,
  age, value`.
- **Indicators:** `tests`, `detections`, `positivity`.
- **Ages:** all ages only.

| File | Holds | Loaded as ([`load_data.R`](../../code/01_main_supporting/load_data.R)) |
|---|---|---|
| [`erviss_detections_sentinel.csv`](erviss_detections_sentinel.csv) | ERVISS `sentinelTestsDetectionsPositivity.csv`: primary-care sentinel swabs | `data$epi$erviss_typing_sentinel` |
| [`erviss_detections_nonsentinel.csv`](erviss_detections_nonsentinel.csv) | ERVISS `nonSentinelTestsDetections.csv`: tests and detections from other sources, mostly hospital laboratories | `data$epi$erviss_typing_nonsentinel` |

**How the ILI+ panel uses them.** The panel takes positivity as detections ÷ tests from the
**sentinel** file. For Croatia, Finland, Iceland, Latvia, Malta and Romania it uses the
**non-sentinel** file instead, because their sentinel testing is absent or too sparse (see
[`stitch_iliplus.R`](../../code/01_main_supporting/stitch_iliplus.R)). Whether that is the right
choice is [open question 4](../../documentation/to_confirm_with_surveillance.md).

## Coverage

Sentinel influenza tests: weeks reported per season, and the median number of tests in a winter week
(weeks 40–20).

| Country | 2021/22 | 2022/23 | 2023/24 | 2024/25 | 2025/26 | Median tests a winter week |
|---|---:|---:|---:|---:|---:|---:|
| Austria | 27 | 21 | 52 | 52 | 33 | 168 |
| Belgium | 41 | 52 | 52 | 52 | 32 | 16 |
| Bulgaria | 51 | 51 | 52 | 50 | 31 | 7 |
| Croatia | 48 | 52 | 52 | 22 | | 0 |
| Cyprus | | | 39 | 52 | 33 | 65 |
| Czechia | 52 | 52 | 52 | 52 | 13 | 42 |
| Denmark | 52 | 52 | 52 | 52 | 33 | 149 |
| Estonia | 49 | 50 | 52 | 52 | 20 | 26 |
| Finland | 52 | 38 | 11 | | | 2 |
| France | 52 | 52 | 52 | 50 | 32 | 115 |
| Germany | 52 | 52 | 52 | 52 | 33 | 149 |
| Greece | 33 | 52 | 52 | 52 | 33 | 54 |
| Hungary | 29 | 33 | 31 | 51 | 33 | 116 |
| Iceland | 52 | 52 | 52 | | | 0 |
| Ireland | 52 | 51 | 52 | 52 | 33 | 99 |
| Italy | 42 | 24 | 24 | 24 | 20 | 250 |
| Latvia | 49 | 46 | 26 | | | 0 |
| Lithuania | 52 | 52 | 51 | 26 | 23 | 14 |
| Luxembourg | 51 | 52 | 52 | 52 | 33 | 62 |
| Malta | | 51 | 45 | | | 0 |
| Netherlands | 52 | 52 | 52 | 52 | 33 | 30 |
| Norway | 52 | 52 | 52 | 52 | 33 | 56 |
| Poland | 51 | 51 | 52 | 51 | 33 | 52 |
| Portugal | 52 | 52 | 51 | 52 | 33 | 19 |
| Romania | 46 | | 32 | 6 | | 2 |
| Slovakia | 52 | 52 | 52 | 50 | 32 | 36 |
| Slovenia | 52 | 52 | 52 | 52 | 33 | 32 |
| Spain | | 52 | 52 | 52 | 33 | 796 |
| Sweden | 49 | 33 | 33 | 33 | 32 | 13 |

**Countries that rely on non-sentinel testing.** Five of the six countries above have little or no
sentinel testing; their median non-sentinel tests in a winter week are:

| Country | Non-sentinel tests |
|---|---:|
| Croatia | 300 |
| Iceland | 216 |
| Latvia | 46 |
| Malta | 542 |
| Romania | 329 |

Finland runs about 2 tests a week in either stream and stops in 2023/24.

## Earlier seasons

**In this repository.** This snapshot starts at 2021-W25. For earlier seasons the project relies on
RespiCompass's ILI+, [`data_respicompass_iliplus.csv`](../data_respicompass_iliplus.csv), in which
positivity is already combined with the ILI rate.

**A lead within ERVISS, not yet used.** ERVISS's own dated snapshots from early to October 2024 carry
sentinel positivity back to 2014-W40: checked for 2024-02-09 to 2024-10-18. The ILI rates are in the
same position, as the [ILI rates folder](../primary_care_ili_rates/README.md#earlier-seasons)
describes.

## Related

- [Primary care ILI rates](../primary_care_ili_rates/README.md): ILI+ is ILI rate × this positivity.
- [SARI test positivity](../sari_flu_positivity/README.md): the same measure among hospitalised SARI
  patients.
- [Open questions 1 and 4](../../documentation/to_confirm_with_surveillance.md): what an absent
  detections row means, and the non-sentinel choice.
