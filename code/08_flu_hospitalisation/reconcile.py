"""How the national series compare with ERVISS for the same country, over the periods both cover.

  DK  SSI laboratory-confirmed admissions   vs  ERVISS hospital admissions (both counts)
  DE  RKI influenza-coded SARI rate         vs  ERVISS SARI+ (both per 100 000)
  ES  SiVIRA SARI+                          vs  ERVISS SARI+ (both per 100 000)
  FR  OSCOUR share after emergency visits   vs  ERVISS ICU admissions (different units: r only)
  SE  monthly J09-J11 inpatient episodes    vs  ERVISS ICU admissions summed by month (both counts)
Austria is not compared: ERVISS has its all-cause SARI rate but no influenza testing. Poland is not in
ERVISS's hospital data.

Usage: python3 -I reconcile.py <out dir>   -> <out dir>/reconciliation_national_vs_erviss.csv
"""
import csv, collections, datetime, statistics, sys

OUT = sys.argv[1]
PAIRS = [('DK', 'Denmark', 'lab_confirmed_admissions', 'hosp_admissions', 'week', True),
         ('DE', 'Germany', 'SARI_influenza_coded_rate', 'sari_plus', 'week', True),
         ('ES', 'Spain', 'SARI_influenza_rate', 'sari_plus', 'week', True),
         ('FR', 'France', 'ED_to_hospital_influenza', 'icu_admissions', 'week', False),
         ('SE', 'Sweden', 'inpatient_episodes_ICD_J09-J11', 'icu_admissions', 'month', True)]

erviss = collections.defaultdict(dict)          # (country, indicator) -> period -> value
for r in csv.DictReader(open(f"{OUT}/erviss_flu_hospital_weekly.csv")):
    if r['age'] == 'total' and r['subtype'] == 'total':
        erviss[(r['country'], r['indicator'])][r['yearweek']] = float(r['value'])
national = collections.defaultdict(dict)
for r in csv.DictReader(open(f"{OUT}/national_flu_hospital_series.csv")):
    national[(r['country'], r['indicator'])][r['period']] = float(r['value'])


def by_month(weekly):   # ISO weeks summed into the calendar month of their Thursday
    m = collections.Counter()
    for yw, v in weekly.items():
        d = datetime.date.fromisocalendar(int(yw[:4]), int(yw.split('W')[1]), 4)
        m[f"{d.year}-{d.month:02d}"] += v
    return m


with open(f"{OUT}/reconciliation_national_vs_erviss.csv", 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['country', 'national_indicator', 'erviss_indicator', 'period_type', 'periods_compared',
                'first_period', 'last_period', 'pearson_r', 'ratio_national_to_erviss'])
    for cc, name, nat, erv, ptype, same_units in PAIRS:
        a = national[(cc, nat)]
        b = by_month(erviss[(name, erv)]) if ptype == 'month' else erviss[(name, erv)]
        common = sorted(set(a) & set(b))
        x, y = [a[p] for p in common], [b[p] for p in common]
        r = statistics.correlation(x, y)
        ratio = round(sum(x) / sum(y), 3) if same_units and sum(y) else ''
        w.writerow([cc, nat, erv, ptype, len(common), common[0], common[-1], round(r, 3), ratio])
        print(f"{cc} {nat:32s} vs {erv:16s} {len(common):4d} {ptype}s  r = {r:.3f}  ratio = {ratio}")
