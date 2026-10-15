"""Standardise the seven machine-readable national influenza-hospitalisation series into one long table.
The inputs are the copies fetch_sources.sh pins (<work_dir>/national/<cc>/...):

  DK  SSI influenza dashboard: new admissions with laboratory-confirmed influenza A+B, weekly
  AT  SARI-Dashboard: SARI admissions with an influenza diagnosis code, weekly (ward + ICU)
  DE  RKI ICOSARI: SARI incidence with ICD-10 J10 coded, per 100 000, weekly (the series ERVISS shows)
  ES  ISCIII SiVIRA: SARI rate x laboratory influenza positivity, per 100 000, weekly (as ERVISS, from 2020/21)
  FR  SpF Odisse (OSCOUR): hospitalisations after an emergency visit for influenza (clinical diagnosis),
      per 100 000 hospitalisations after an emergency visit, weekly
  SE  Socialstyrelsen: inpatient care episodes with main diagnosis J09-J11, MONTHLY
  PL  NIZP-PZH Meldunki: persons referred to hospital for influenza or suspected influenza, quarter-months

Usage: python3 -I build_national.py <work_dir> <out dir>   -> <out dir>/national_flu_hospital_series.csv
"""
import csv, collections, datetime, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # xlsx_read.py sits beside this file
from xlsx_read import read

W, OUT = sys.argv[1:3]
N = f"{W}/national"
out = []   # (country, source, indicator, unit, period_type, period, mid_date, season, value)


def season_of(d):
    return f"{d.year}/{d.year + 1}" if d.month >= 8 else f"{d.year - 1}/{d.year}"


def add_week(c, src, ind, unit, yw, v):
    y, w = int(yw[:4]), int(yw.split("W")[1])
    d = datetime.date.fromisocalendar(y, w, 4)
    out.append((c, src, ind, unit, "week", yw, d.isoformat(), season_of(d), v))


# DK: region 'Alle', sex 'Alle', INFLA + INFLB summed over the age groups
dk = collections.Counter()
for r in csv.DictReader(open(f"{N}/dk/docs/02_influenza_epikurve_season_region_uge_agegrp.csv", encoding="utf-8-sig"), delimiter=";"):
    if r["Region"] == "Alle" and r["Køn"] == "Alle" and r["Aldersgruppe"] != "Alle":
        dk[r["Uge"].replace("-U", "-W")] += float(r["Antal nye indlæggelser"] or 0)
for yw, v in dk.items():
    add_week("DK", "SSI influenza dashboard", "lab_confirmed_admissions", "count", yw, v)

# AT: SARI admissions coded influenza, all regions, sexes, ages and wards
at = collections.Counter()
for r in json.load(open(f"{N}/at/sari/patient.json"))["data"]:
    m = re.match(r"(\d+)\. KW (\d{4})", r["KW"])
    at[f"{m.group(2)}-W{int(m.group(1)):02d}"] += r["INFLUENZA"]
for yw, v in at.items():
    add_week("AT", "SARI-Dashboard (Sozialversicherung)", "SARI_influenza_coded_admissions", "count", yw, v)

# DE: influenza-SARI incidence, all ages
for r in csv.DictReader(open(f"{N}/de/SARI-Hospitalisierungsinzidenz.tsv"), delimiter="\t"):
    if r["SARI"] == "Influenza" and r["Altersgruppe"] == "00+":
        add_week("DE", "RKI ICOSARI", "SARI_influenza_coded_rate", "per 100 000", r["Kalenderwoche"], float(r["Hospitalisierungsinzidenz"]))

# ES: influenza-attributable SARI rate (column A = rate, C = year-week)
for r in read(f"{N}/es/tablas_informe_SiVIRA_202551.xlsx", "tasa_gripe_temporadas_nacional")[1:]:
    try:
        add_week("ES", "ISCIII SiVIRA", "SARI_influenza_rate", "per 100 000", f"{r['C'][:4]}-W{r['C'][4:]}", float(r["A"]))
    except (KeyError, TypeError, ValueError):
        pass   # header and note rows

# FR: share of hospitalisations after an emergency visit that are for influenza, all ages
FR_COL = "Taux d'hospitalisations après passages aux urgences pour grippe"
for r in csv.DictReader(open(f"{N}/fr/data/raw/grippe-passages-aux-urgences-et-actes-sos-medecins-france.csv", encoding="utf-8-sig")):
    if r["Classe d'âge"] == "Tous âges" and r[FR_COL]:
        y, w = r["Semaine"].split("-S")
        add_week("FR", "SpF Odissé (OSCOUR)", "ED_to_hospital_influenza", "per 100 000 hospitalisations after an emergency visit",
                 f"{y}-W{int(w):02d}", float(r[FR_COL]))

# SE: monthly J09-J11 main diagnosis, both sexes, all ages
MON = {m: i + 1 for i, m in enumerate("Januari Februari Mars April Maj Juni Juli Augusti September Oktober November December".split())}
se = collections.Counter()
for r in csv.DictReader(open(f"{N}/se/postdata/2024-11-25-vag/data/infl_se_par_085.tsv"), delimiter="\t"):
    se[(int(r["Year"]), MON[r["Month"]])] += sum(int(r[k] or 0) for k in ("J09", "J10", "J11"))
for (y, m), v in se.items():
    d = datetime.date(y, m, 15)
    out.append(("SE", "Socialstyrelsen inpatient register", "inpatient_episodes_ICD_J09-J11", "count", "month",
                f"{y}-{m:02d}", d.isoformat(), season_of(d), v))

# PL: quarter-month periods (days 1-7, 8-15, 16-22, 23-end), labelled A-D
pl = collections.Counter(); pdays = {}
for f in ("g_szpital_2010_2022.csv", "g_szpital_2023.csv"):
    for r in csv.DictReader(open(f"{N}/pl/{f}")):
        mo = int(re.match(r"(\d+)", r["numer"]).group(1)); k = (int(r["rok"]), mo, r["numer"][-1])
        pl[k] += float(r["o_l"] or 0); pdays[k] = r["dni"]
for (y, mo, q), v in pl.items():
    a, b = pdays[(y, mo, q)].split("-")
    d = datetime.date(y, mo, (int(a) + int(b)) // 2)
    out.append(("PL", "NIZP-PZH Meldunki", "hospital_referrals_influenza_or_suspected", "count", "quarter-month",
                f"{y}-{mo:02d}{q}", d.isoformat(), season_of(d), v))

os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/national_flu_hospital_series.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["country", "source", "indicator", "unit", "period_type", "period", "mid_date", "season", "value"])
    w.writerows(sorted(out))
cnt = collections.Counter((o[0], o[2]) for o in out)
for k in sorted(cnt):
    ps = sorted(o[6] for o in out if (o[0], o[2]) == k)
    print(f"{k[0]} {k[1]:42s} {cnt[k]:5d} periods  {ps[0]} -> {ps[-1]}")
