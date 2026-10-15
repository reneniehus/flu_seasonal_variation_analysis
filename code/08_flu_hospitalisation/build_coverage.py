"""Which influenza-hospitalisation indicator is available for each country and season, and which is best.

Every series in the two weekly tables is put in one indicator class; the letter is the preference order.
  A  laboratory-confirmed influenza admissions
  B  SARI rate x laboratory influenza positivity of SARI samples (SARI+)
  C  SARI with an influenza diagnosis code (DE: ICD-10 J10; AT: diagnosis-coded)
  D  hospitalisations after an emergency visit for influenza (clinical diagnosis; FR)
  E  inpatient care episodes with main diagnosis J09-J11, monthly (SE)
  F  ICU or severe cases only
  G  hospital referrals for influenza or suspected influenza (PL)
ERVISS 'hospital admissions' are counted as ICU (F) where they are the ICU series: in Finland,
Netherlands, Portugal and Sweden the two are identical in 85-100% of weeks; in France ICU is 1.00 of
admissions over the seasons and in Italy 0.94. Germany's ERVISS SARI+ is the J10-coded share of SARI
(RKI), not a laboratory positivity, so it is C.

The best class for a country-season is the first, by letter, with at least 20 periods (5 for the
monthly class E). If no class reaches that, it is the class with the most periods.

Usage: python3 -I build_coverage.py <out dir>   (reads and writes there)
"""
import csv, collections, datetime, sys

OUT = sys.argv[1]
ISO = {'Austria': 'AT', 'Belgium': 'BE', 'Bulgaria': 'BG', 'Croatia': 'HR', 'Cyprus': 'CY', 'Czechia': 'CZ',
       'Denmark': 'DK', 'Estonia': 'EE', 'Finland': 'FI', 'France': 'FR', 'Germany': 'DE', 'Greece': 'GR',
       'Hungary': 'HU', 'Iceland': 'IS', 'Ireland': 'IE', 'Italy': 'IT', 'Latvia': 'LV', 'Lithuania': 'LT',
       'Luxembourg': 'LU', 'Malta': 'MT', 'Netherlands': 'NL', 'Norway': 'NO', 'Poland': 'PL', 'Portugal': 'PT',
       'Romania': 'RO', 'Slovakia': 'SK', 'Slovenia': 'SI', 'Spain': 'ES', 'Sweden': 'SE'}
CLASS = {'A': 'A lab-confirmed admissions', 'B': 'B SARI x lab influenza positivity', 'C': 'C SARI coded as influenza',
         'D': 'D emergency visits leading to admission', 'E': 'E inpatient diagnoses J09-J11 (monthly)',
         'F': 'F ICU / severe cases only', 'G': 'G hospital referrals, suspected influenza'}
ICU_LIKE = {'Finland', 'France', 'Italy', 'Netherlands', 'Portugal', 'Sweden'}
NATIONAL = {'lab_confirmed_admissions': 'A', 'SARI_influenza_rate': 'B', 'SARI_influenza_coded_rate': 'C',
            'SARI_influenza_coded_admissions': 'C', 'ED_to_hospital_influenza': 'D',
            'inpatient_episodes_ICD_J09-J11': 'E', 'hospital_referrals_influenza_or_suspected': 'G'}


def season(d):   # seasons run August to July, as in the rest of the repository
    return f"{d.year}/{d.year + 1}" if d.month >= 8 else f"{d.year - 1}/{d.year}"


def erviss_class(indicator, country):
    if indicator == 'hosp_admissions':
        return 'F' if country in ICU_LIKE else 'A'
    if indicator == 'icu_admissions':
        return 'F'
    if indicator == 'sari_plus':
        return 'C' if country == 'Germany' else 'B'
    return None


cells = collections.defaultdict(lambda: collections.defaultdict(set))   # (cc, season) -> class -> periods
sources = collections.defaultdict(set)                                  # (cc, season, class) -> sources
for r in csv.DictReader(open(f"{OUT}/erviss_flu_hospital_weekly.csv")):
    if r['age'] != 'total' or r['subtype'] != 'total' or r['country'] not in ISO:
        continue
    cl = erviss_class(r['indicator'], r['country'])
    if cl is None:
        continue
    y, w = int(r['yearweek'][:4]), int(r['yearweek'].split('W')[1])
    key = (ISO[r['country']], season(datetime.date.fromisocalendar(y, w, 4)))
    cells[key][cl].add(r['yearweek'])
    sources[key + (cl,)].add('ERVISS')
for r in csv.DictReader(open(f"{OUT}/national_flu_hospital_series.csv")):
    cl = NATIONAL[r['indicator']]
    key = (r['country'], season(datetime.date.fromisoformat(r['mid_date'])))
    cells[key][cl].add(r['period'])
    sources[key + (cl,)].add(r['source'])

with open(f"{OUT}/coverage_by_country_season.csv", 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['country', 'season', 'best_class', 'periods', 'best_sources', 'all_classes'])
    for (cc, s), d in sorted(cells.items()):
        n = {cl: len(p) for cl, p in d.items()}
        full = [cl for cl in sorted(n) if n[cl] >= (5 if cl == 'E' else 20)]
        best = full[0] if full else max(sorted(n), key=lambda cl: n[cl])
        w.writerow([cc, s, CLASS[best], n[best], '; '.join(sorted(sources[(cc, s, best)])),
                    ' | '.join(f"{CLASS[cl]} ({n[cl]})" for cl in sorted(n))])
print('country-seasons', len(cells), '| countries', len({k[0] for k in cells}))
