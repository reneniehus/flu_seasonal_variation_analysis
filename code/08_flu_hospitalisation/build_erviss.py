"""Clean weekly influenza-hospitalisation series from ERVISS, merged across the current files and the
dated snapshots, so each series reaches as far back as ERVISS has ever published it.

Merge rule: for every key (indicator, country, age, subtype, week) the value comes from the NEWEST file
that contains it: the current file first, then the snapshots from newest to oldest. Older snapshots
therefore only fill weeks that newer files have dropped. (The files are a rolling window; the 2024
severity snapshots still carry 2014/15-2019/20.)

Indicators (influenza only):
  sari_rate                SARI syndromic rate per 100 000 (sentinel hospitals), by age
  sari_flu_tests, sari_flu_detections, sari_flu_positivity
                           SARI virological, all ages; detections also by (sub)type
  sari_plus                sari_rate x influenza positivity / 100. Positivity is all-age, so it is
                           assumed equal across ages (owner decision, 2026-10).
  hosp_admissions, icu_admissions, icu_inpatients, hospital_inpatients, deaths
                           non-sentinel severity counts (laboratory-confirmed influenza), by age

Usage: python3 -I build_erviss.py <erviss repo>/data <out dir>
"""
import csv, glob, os, sys, collections

DATA, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)


def files(stem):
    snaps = sorted(glob.glob(f"{DATA}/snapshots/*_{stem}.csv"), reverse=True)   # newest first
    return [(f"{DATA}/{stem}.csv", "current")] + [(p, os.path.basename(p)[:10]) for p in snaps]


def season(yw):
    y, w = int(yw[:4]), int(yw.split("W")[1])
    return f"{y}/{y + 1}" if w >= 40 else f"{y - 1}/{y}"


def num(v):
    return None if v in ("", "NA", None) else float(v)


merged = {}                                    # key -> (value, source file)
revisions = collections.defaultdict(list)      # indicator -> relative differences between files


def take(key, val, src, indicator):
    if val is None:
        return
    if key in merged:
        v0 = merged[key][0]
        if v0 != val and max(abs(v0), abs(val)) > 0:
            revisions[indicator].append(abs(v0 - val) / max(abs(v0), abs(val)))
        else:
            revisions[indicator].append(0.0)
        return
    merged[key] = (val, src)


SEV = {"hospitaladmissions": "hosp_admissions", "ICUadmissions": "icu_admissions",
       "ICUinpatients": "icu_inpatients", "hospitalinpatients": "hospital_inpatients", "deaths": "deaths"}
for path, src in files("nonSentinelSeverity"):
    for r in csv.DictReader(open(path)):
        if r["pathogen"] != "Influenza" or r["indicator"] not in SEV:
            continue
        ind = SEV[r["indicator"]]
        take((ind, r["countryname"], r["age"], "total", r["yearweek"]), num(r["value"]), src, ind)
for path, src in files("SARIRates"):
    for r in csv.DictReader(open(path)):
        take(("sari_rate", r["countryname"], r["age"], "total", r["yearweek"]), num(r["value"]), src, "sari_rate")
for path, src in files("SARITestsDetectionsPositivity"):
    for r in csv.DictReader(open(path)):
        if r["pathogen"] != "Influenza":
            continue
        sub = r["pathogensubtype"] if r["pathogensubtype"] not in ("NA", "") else r["pathogentype"]
        ind = "sari_flu_" + r["indicator"]
        take((ind, r["countryname"], r["age"], sub, r["yearweek"]), num(r["value"]), src, ind)

# SARI+ = SARI rate (each age group) x all-age influenza positivity of SARI samples
pos = {}
for (ind, cn, age, sub, yw), (v, src) in merged.items():
    if sub == "total" and age == "total" and ind == "sari_flu_tests":
        det = merged.get(("sari_flu_detections", cn, "total", "total", yw))
        p = merged.get(("sari_flu_positivity", cn, "total", "total", yw))
        if p is not None:
            pos[(cn, yw)] = (p[0], v)
        elif det is not None and v > 0:
            pos[(cn, yw)] = (100 * det[0] / v, v)
for (ind, cn, age, sub, yw), (v, src) in list(merged.items()):
    if ind == "sari_rate" and (cn, yw) in pos:
        p, tests = pos[(cn, yw)]
        merged[("sari_plus", cn, age, "total", yw)] = (v * p / 100, src)

with open(f"{OUT}/erviss_flu_hospital_weekly.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["indicator", "country", "season", "yearweek", "age", "subtype", "value", "unit", "source_file", "flu_tests_that_week"])
    unit = {"sari_rate": "per 100 000", "sari_plus": "per 100 000", "sari_flu_positivity": "%"}
    for (ind, cn, age, sub, yw), (v, src) in sorted(merged.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][4], kv[0][2], kv[0][3])):
        tests = pos.get((cn, yw), (None, None))[1] if ind in ("sari_plus", "sari_flu_positivity") else None
        w.writerow([ind, cn, season(yw), yw, age, sub, round(v, 6), unit.get(ind, "count"), src, "" if tests is None else int(tests)])

# coverage: weeks with data per country x season, all ages, for every indicator
cov = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
for (ind, cn, age, sub, yw), (v, src) in merged.items():
    if age == "total" and sub == "total":
        cov[ind][cn][season(yw)] += 1
with open(f"{OUT}/erviss_flu_hospital_coverage.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["indicator", "country", "season", "weeks"])
    for ind in sorted(cov):
        for cn in sorted(cov[ind]):
            for s in sorted(cov[ind][cn]):
                w.writerow([ind, cn, s, cov[ind][cn][s]])
for ind in sorted(revisions):
    d = revisions[ind]
    big = sum(x > 0.05 for x in d)
    print(f"revision check {ind:22s}: {len(d):6d} overlapping values across files, {big:5d} differ by >5% ({100 * big / max(len(d), 1):.1f}%)")
print("keys written:", len(merged))
