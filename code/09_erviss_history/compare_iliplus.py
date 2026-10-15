"""Does RespiCompass's influenza ILI+ equal ERVISS's own ILI rate x positivity?

The ILI+ panel takes 2014/15-2023/24 from RespiCompass (data/data_respicompass_iliplus.csv). The
merged ERVISS releases now reach back to 2014-W40 as well, so the two can be compared week by week.
ERVISS ILI+ = ILI consultation rate (all ages) x influenza detections / tests (all ages). Each
RespiCompass row says which tests it used ('sentinel' or 'non_sentinel'), and the same stream is used
here, so the comparison is like for like.

Writes <data>/primary_care_ili_rates/iliplus_erviss_vs_respicompass.csv, per country:
  weeks compared, first and last week, Pearson r,
  median ratio RespiCompass / ERVISS (weeks where both are above zero),
  share of weeks agreeing within 1%.

Usage: python3 -I compare_iliplus.py <repository data dir>
"""
import collections, csv, statistics, sys

DATA = sys.argv[1]
ili = {}
for r in csv.DictReader(open(f"{DATA}/primary_care_ili_rates/erviss_iliari_merged.csv")):
    if r["indicator"] == "ILIconsultationrate" and r["age"] == "total":
        ili[(r["countryname"], r["yearweek"])] = float(r["value"])
tests = {}
for stream, f in (("sentinel", "erviss_detections_sentinel_merged.csv"), ("non_sentinel", "erviss_detections_nonsentinel_merged.csv")):
    t, d = {}, {}
    for r in csv.DictReader(open(f"{DATA}/primary_care_flu_positivity/{f}")):
        if r["pathogentype"] == "Influenza" and r["age"] == "total" and r["indicator"] in ("tests", "detections"):
            (t if r["indicator"] == "tests" else d)[(r["countryname"], r["yearweek"])] = float(r["value"])
    tests[stream] = {k: d[k] / t[k] for k in t if t[k] > 0 and k in d}

pairs = collections.defaultdict(list)   # country -> (yearweek, respicompass, erviss)
for r in csv.DictReader(open(f"{DATA}/data_respicompass_iliplus.csv")):
    if r["age"] != "total" or r["value"] in ("", "NA"):
        continue
    k = (r["location_name"], r["yearweek"])
    pos = tests[r["source"]].get(k)
    if pos is not None and k in ili:
        pairs[r["location_name"]].append((r["yearweek"], float(r["value"]), ili[k] * pos))

with open(f"{DATA}/primary_care_ili_rates/iliplus_erviss_vs_respicompass.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["country", "weeks_compared", "first_week", "last_week", "pearson_r",
                "median_ratio_respicompass_to_erviss", "share_within_1pct"])
    for c in sorted(pairs):
        p = sorted(pairs[c])
        x, y = [a for _, a, _ in p], [b for _, _, b in p]
        r = statistics.correlation(x, y) if len(p) > 2 and len(set(x)) > 1 and len(set(y)) > 1 else float("nan")
        ratios = [a / b for _, a, b in p if a > 0 and b > 0]
        within = sum(abs(a - b) <= 0.01 * max(abs(a), abs(b)) for _, a, b in p) / len(p)
        med = statistics.median(ratios) if ratios else float("nan")
        w.writerow([c, len(p), p[0][0], p[-1][0], round(r, 3), round(med, 3), round(within, 3)])
        print(f"{c:12s} {len(p):4d} weeks {p[0][0]}..{p[-1][0]}  r = {r:.3f}  ratio = {med:.3f}  within 1%: {within:.0%}")
