"""ERVISS's primary-care files merged across its releases, so each series reaches back as far as ERVISS
has ever published it: 2014-W40.

Merge rule: for every key (all columns but `value`), the value comes from the NEWEST release that has
it: the current file, then the releases of 2026-05-22, 2025-06-20 and 2024-10-18 (fetch_releases.sh).
A week therefore carries its last-published value. This is the same rule as the hospital indicators
(code/08_flu_hospitalisation/build_erviss.py).

Outputs, in the repository's data folders:
  primary_care_ili_rates/erviss_iliari_merged.csv                        ILI and ARI rates, all ages
  primary_care_flu_positivity/erviss_detections_sentinel_merged.csv      influenza rows only
  primary_care_flu_positivity/erviss_detections_nonsentinel_merged.csv   influenza rows only
Each keeps ERVISS's columns and adds two:
  release        the release the value came from ('current' or its date)
  dropped_later  TRUE when a newer release covers that week but no longer carries the value. Such values
                 are kept, so the user can choose; whether ECDC withdrew them on purpose is not known.
Both folders also get coverage_by_country_season.csv: weeks with data (for tests, weeks with at least one
test); seasons run August to July, as everywhere in the repository.

Printed checks, per file:
  revised  share of values present in two releases that changed by more than 5%
  filled   keys taken from an older release for a week that a newer release also covers, i.e. values
           a newer release no longer carries

Usage: python3 -I merge_releases.py <erviss repo>/data <repository data dir>
"""
import collections, csv, datetime, glob, os, sys

SRC, DATA = sys.argv[1:3]
FILES = [  # ERVISS stem, output, pathogen kept (None = all rows)
    ("ILIARIRates", "primary_care_ili_rates/erviss_iliari_merged.csv", None),
    ("sentinelTestsDetectionsPositivity", "primary_care_flu_positivity/erviss_detections_sentinel_merged.csv", "Influenza"),
    ("nonSentinelTestsDetections", "primary_care_flu_positivity/erviss_detections_nonsentinel_merged.csv", "Influenza"),
]


def season(yw):   # August to July
    d = datetime.date.fromisocalendar(int(yw[:4]), int(yw.split("W")[1]), 4)
    return f"{d.year}/{d.year + 1}" if d.month >= 8 else f"{d.year - 1}/{d.year}"


def releases(stem):   # newest first
    snaps = sorted(glob.glob(f"{SRC}/snapshots/*_{stem}.csv"), reverse=True)
    return [(f"{SRC}/{stem}.csv", "current")] + [(p, os.path.basename(p)[:10]) for p in snaps]


merged_by_file = {}
for stem, out, pathogen in FILES:
    merged, columns, windows = {}, None, []
    revised = changed = filled = 0
    filled_by_country = collections.Counter()
    for path, rel in releases(stem):
        rows = [r for r in csv.DictReader(open(path)) if pathogen is None or r["pathogen"] == pathogen]
        if columns is None:
            columns = list(rows[0].keys())
        key_cols = [c for c in columns if c != "value"]
        weeks = sorted({r["yearweek"] for r in rows})
        for r in rows:
            if r["value"] in ("", "NA"):
                continue
            key = tuple(r[c] for c in key_cols)
            if key in merged:
                v0, v1 = float(merged[key][0]), float(r["value"])
                revised += 1
                changed += abs(v0 - v1) > 0.05 * max(abs(v0), abs(v1))
                continue
            dropped = any(start <= r["yearweek"] <= end for start, end in windows)
            if dropped:
                filled += 1
                filled_by_country[r["countryname"]] += 1
            merged[key] = (r["value"], rel, dropped)
        windows.append((weeks[0], weeks[-1]))
        print(f"{stem:35s} {rel:10s} {weeks[0]} to {weeks[-1]}  {len(rows):7d} rows")
    os.makedirs(os.path.dirname(f"{DATA}/{out}"), exist_ok=True)
    with open(f"{DATA}/{out}", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns + ["release", "dropped_later"])
        yi, ci = key_cols.index("yearweek"), key_cols.index("countryname")
        for key in sorted(merged, key=lambda k: (k[ci], k[yi], k)):
            row = dict(zip(key_cols, key))
            value, rel, dropped = merged[key]
            w.writerow([value if c == "value" else row[c] for c in columns] + [rel, "TRUE" if dropped else "FALSE"])
    print(f"  -> {out}: {len(merged)} values; revised: {changed} of {revised} shared values changed by >5% "
          f"({100 * changed / max(revised, 1):.1f}%); filled from older releases: {filled} {dict(filled_by_country.most_common(5))}")
    merged_by_file[stem] = (merged, key_cols)

# coverage: weeks with data per country and season (all ages; for the tests files, weeks with at least one
# influenza test, since a week of zero tests has no positivity)
cov = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
for stem, (merged, key_cols) in merged_by_file.items():
    for key in merged:
        k = dict(zip(key_cols, key))
        if k["age"] != "total":
            continue
        if stem == "ILIARIRates":
            col = "ili_weeks" if k["indicator"] == "ILIconsultationrate" else "ari_weeks"
        elif k["pathogentype"] == "Influenza" and k["indicator"] == "tests" and float(merged[key][0]) > 0:
            col = "sentinel_test_weeks" if stem.startswith("sentinel") else "nonsentinel_test_weeks"
        else:
            continue
        cov["ili" if stem == "ILIARIRates" else "pos"][(k["countryname"], season(k["yearweek"]))][col] += 1
for name, folder, cols in (("ili", "primary_care_ili_rates", ["ili_weeks", "ari_weeks"]),
                           ("pos", "primary_care_flu_positivity", ["sentinel_test_weeks", "nonsentinel_test_weeks"])):
    with open(f"{DATA}/{folder}/coverage_by_country_season.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "season"] + cols)
        for (cn, s), c in sorted(cov[name].items()):
            w.writerow([cn, s] + [c[x] for x in cols])
print("coverage tables written")
