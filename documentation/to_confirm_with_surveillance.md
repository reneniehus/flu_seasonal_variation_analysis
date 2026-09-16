# To confirm with surveillance colleagues

Questions about **what the source data mean** that we have answered by inference and cannot answer from
the data alone. Each entry states the question, what we currently assume, the evidence behind that
assumption, and **what would change if the answer is different** — so an answer can be acted on without
re-deriving the analysis.

These are not code defects. They are places where the model rests on a reading of a publication format
or a surveillance practice, and where being wrong would change a result rather than break a run.

Status key: **OPEN** (asked, no answer) · **ANSWERED** (record the answer and the date) · **CLOSED**
(answered and acted on).

---

## 1. Does an absent `detections` row mean zero detections, or an unpublished count? — **OPEN**

**The question.** In the ERVISS typing files a week with no influenza detections is encoded two
different ways. **3151 weeks state `detections = 0` explicitly** (2083 sentinel, 1068 non-sentinel).
**709 weeks report `tests > 0` with the detections row simply absent** (271 sentinel, 438
non-sentinel). No published `positivity` value rescues any of the 709. Which is it?

**Why it matters to us.** We re-derive positivity as detections/tests, so an explicit zero becomes an
observed zero while an absent row becomes missing — and a missing ILI+ week is dropped from the panel.
Because the season grid ends at the last finite week, a *trailing* run of them shortens the season
rather than leaving holes. Across the whole 25-country panel this affects **301 weeks in 34
country-seasons of 12 countries** (AT, BG, CZ, HU, IS, IT, LT, LV, MT, PL, RO, SK).

**What we currently assume.** That an absent row is *not* reliably a zero, so a country-season is
excluded when it contains weeks that cannot plausibly have been zero. In the fitted design that
excludes **CZ 2024/2025** and keeps everything else, giving 12 countries / 85 country-seasons / 182
parameters.

**The evidence for that reading** (analysed 2026-09-16; `erviss_encoding_ambiguous()` in
`code/01_main_supporting/stitch_iliplus.R`):

- **Omission is not the format's way of writing zero.** 18 of 30 countries have **zero** absent rows
  and write explicit zeros throughout. Only IT and MT omit without ever writing an explicit zero. Ten
  countries do **both** — CZ has 79 explicit zeros *and* 30 absent rows; PL has 136 and 1. If omission
  meant zero, the countries writing thousands of explicit zeros would not also omit.
- **Absent weeks are not preferentially off-season**: median season week 19 with 32% at week 40 or
  later, against 24 and 37% for explicit zeros.
- **Per-week plausibility.** For each affected week we compute P(0 detections | that week's test count,
  the local positivity of published weeks within ±3 weeks) under a binomial. Zero is plausible
  (P ≥ 0.05) for **82%** of the 177 weeks that have usable neighbours, against **96%** for genuine
  explicit zeros. So most absent rows do look like quiet weeks — but a minority demonstrably cannot be,
  and 10 are essentially impossible as zeros.
- **The two cases in the fitted design differ in kind**, which is why the exclusion is decided per week
  rather than per season:
  - **CZ 2024/2025** — week 34 has 56 tests against 7.3% local positivity, so about 4 detections were
    expected (P(0) = 0.014); weeks 37–39 are plausible; **weeks 40–52 have no published neighbour at
    all**, because CZ's detections feed went dark on 2025-03-26 and resumed only the following season.
    14 of its 17 affected weeks cannot plausibly have been zero. Excluded.
  - **PL 2024/2025** — one affected week, season week 3, 13 tests, neighbours showing **0 detections
    over 114 tests** (P(0) = 1.000), flanked by explicit zeros. Almost certainly a genuine zero. Kept.

**What changes with the answer.**

- **If an absent row means zero detections:** the fix is upstream in `code/01_main_supporting/gen_model_input.R`
  — read the absent row as 0 rather than NA — and **this exclusion should be reverted**, recovering
  110 panel weeks and CZ 2024/2025. Note this would also restore CZ's ERVISS baseline slot, which is
  currently gone because 2024/2025 was CZ's only ERVISS season.
- **If it means the count is genuinely unpublished:** the exclusion stays, and a *second* defect still
  needs fixing — a trailing run of missing weeks shortens the season grid instead of leaving holes, so
  a country-season can be fitted on a window that contains no off-season with nothing reporting it.
- **If the convention differs by country** (IT and MT never write an explicit zero, so omission may
  genuinely be their convention while being a lapse elsewhere): the rule should be per country, and
  that matters for any future expansion of the design — Iceland and Malta both otherwise qualify for
  inclusion.

**Where this is implemented.** `erviss_encoding_ambiguous()` derives the affected weeks and their
plausibility from `models_in`; `jm_build_data(exclude_ambiguous_positivity = FALSE)` fits everything
anyway; `ambiguous_min_unexplained` raises the tolerance. Both designs are pinned by the test suite.

---

## 2. Are the two ILI+ sources' off-season floors comparable at all? — **OPEN**

**The question.** We fit one off-season baseline per data source per country (`b_c,src`), and the
recorded reason was that RespiCompass ILI+ is exactly zero in weeks without detections while the ERVISS
reconstruction sits on a small positive floor. **That is contradicted by the data**: over the 357
country-weeks of the 2023/2024 overlap where both sources give a finite value, RespiCompass is zero
with ERVISS positive **0 times**, and the reverse **0 times**; they are zero together 112 times. Where
both are observable they agree about zeros week for week.

**What we now assume.** That the fitted gap between the two baselines (up to 392× for BE) reflects
**which weeks each source covers** — pre-COVID RespiCompass series stop around week 36–42 while
ERVISS-era series run to week 52–53, so one baseline is fitted largely to the wave and the other
largely to the off-season. We therefore read `b` as a coverage-era nuisance parameter, not as a
property of a surveillance system.

**What to confirm.** Is there any real difference in how the two reconstructions encode a quiet week
that we should be modelling, or is our coverage explanation right? And: is the per-country alignment
factor (median RespiCompass/ERVISS over the 2023/2024 overlap) the right way to put them on one scale?

**What changes with the answer.** If the sources genuinely differ in their floor, `b` per source is
substantive and should be interpreted; if not, it stays a nuisance parameter and should never be
reported as a surveillance finding.

---

## 3. Do the pre-COVID RespiCompass series really end around week 36–42? — **OPEN**

**The question.** Pre-COVID RespiCompass country-seasons stop at season week 36–42 while ERVISS-era
seasons run to 52–53. Is that the real extent of reporting, or is our copy of the data truncated?

**Why it matters.** It is the mechanism behind question 2, and it is why some countries' off-season
baselines are estimated from almost no off-season. It also drives the uneven season support the model
now reports (`jm_summary_season`).

---

## 4. Is non-sentinel positivity the right choice for HR, IS, MT, RO, LV and FI? — **OPEN**

**The question.** For most countries our ERVISS ILI+ uses **sentinel** positivity; for these six it
uses **non-sentinel**, because sentinel data are absent or too sparse. Is that the right call for each,
and does it make their ILI+ comparable with the sentinel-based countries?

**Why it matters.** Croatia is in the fitted design. The non-sentinel stream is also the one with the
*more* ambiguous detections rows (438 against 271), so question 1 bites harder there — and Iceland and
Malta, both non-sentinel, are the countries that would otherwise qualify for inclusion.

---

## 5. Is the 2025/2026 season complete? — **OPEN**

**The question.** In our data 2025/2026 runs to 2026-05-13 and only 7 of the 12 design countries report
it at all, on grids of 34–41 weeks against 52–53 for other seasons. Is that season finished as far as
ERVISS is concerned, or still being reported?

**Why it matters.** 2025/2026 carries the **highest** fitted transmissibility of the eight seasons
(R0 1.71) on roughly half the weekly data of the richest season, and the average-one constraint on the
season visibility deviations treats it as one of eight equals. If the season is incomplete rather than
short, it should arguably be excluded or down-weighted rather than compared like for like.

---

## 6. Is the ×1000 rescale for CY, LU and MT correct? — **OPEN**

**The question.** These three report ILI per 100 consultations rather than per 100 000 population, and
we multiply by 1000 to put them on the common basis. Is that the right conversion, and is the resulting
series comparable with a population-denominator rate?

**Why it matters.** None of the three is in the current design, but Malta is an expansion candidate,
and a wrong basis would put its reporting level and age offsets on a different scale from everyone
else's without anything flagging it.
