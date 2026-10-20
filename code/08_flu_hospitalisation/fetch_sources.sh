#!/usr/bin/env bash
# Fetch the exact inputs of the influenza-hospitalisation build, pinned to the commits used in
# October 2026. Blob-less sparse clones, so only the files the build reads are downloaded.
#
# Usage: bash code/08_flu_hospitalisation/fetch_sources.sh <work_dir>
#
# Only github.com is needed. The official hosts (SSI, sozialversicherung.at, ISCIII, Odisse,
# Socialstyrelsen, PZH) were unreachable from the build environment, so the national series are read
# from public GitHub copies; each line names the official source its copy reproduces.
set -euo pipefail
W=${1:?usage: fetch_sources.sh <work_dir>}
mkdir -p "$W"

fetch() {   # fetch <dir> <repo url> <commit> <path>...
  local dir=$W/$1 url=$2 commit=$3; shift 3
  [ -d "$dir/.git" ] || git clone -q --filter=blob:none --no-checkout "$url" "$dir"
  git -C "$dir" sparse-checkout set --no-cone "$@"
  git -C "$dir" checkout -q "$commit"
}

# ERVISS (ECDC): the current files plus dated snapshots, one about every eight weeks. The files are a
# rolling window, so older weeks survive only in older snapshots (build_erviss.py merges them).
SARI_SNAPS="2024-11-08 2024-11-15 2025-01-10 2025-03-07 2025-05-09 2025-07-11 2025-09-05 2025-10-31
            2026-01-09 2026-03-06 2026-05-08 2026-09-25"
SEV_SNAPS="2023-11-24 2023-12-01 2024-01-26 2024-03-22 2024-05-31 2024-07-26 2024-09-20 2024-11-22
           2025-01-24 2025-03-21 2025-05-23 2025-07-25 2025-09-19 2025-11-14 2026-01-23 2026-03-20
           2026-08-07 2026-10-09"
paths=(/data/nonSentinelSeverity.csv /data/SARIRates.csv /data/SARITestsDetectionsPositivity.csv)
for d in $SARI_SNAPS; do paths+=("/data/snapshots/${d}_SARIRates.csv" "/data/snapshots/${d}_SARITestsDetectionsPositivity.csv"); done
for d in $SEV_SNAPS; do paths+=("/data/snapshots/${d}_nonSentinelSeverity.csv"); done
fetch erviss https://github.com/EU-ECDC/Respiratory_viruses_weekly_data c41550037d19ede7988bfa23ea4261abfc28b7cc "${paths[@]}"

# DK: Statens Serum Institut influenza dashboard data (zip 'infl'), mirrored weekly
fetch national/dk https://github.com/steenhulthin/infectious-diseases-data 95cfd34f13439c29f00eca3a448a502fe17d43d2 \
  /docs/02_influenza_epikurve_season_region_uge_agegrp.csv
# AT: SARI-Dashboard open data, opendata-files.sozialversicherung.at/sari/SARI_Wohnregion_Patient_v202602.csv, as JSON
fetch national/at https://github.com/pasrom/wastewater 37f7cfd1182a3749bf070cbc3f04a216e2e239ce /sari/patient.json
# DE: Robert Koch-Institut, ICOSARI (official repository)
fetch national/de https://github.com/robert-koch-institut/SARI-Hospitalisierungsinzidenz aeea12549e46df1f51a92efccc4afa04f5673342 \
  /SARI-Hospitalisierungsinzidenz.tsv
# ES: ISCIII SiVIRA weekly report tables, docsivira.isciii.es/excel/tablas_informe_SiVIRA_202551.xlsx
fetch national/es https://github.com/ccarballo50/sivira-trends-lag f662e481d736d72dd8bce7dab005b0af5e68efba \
  /tablas_informe_SiVIRA_202551.xlsx
# FR: Sante publique France, Odisse dataset grippe-passages-aux-urgences-et-actes-sos-medecins-france
fetch national/fr https://github.com/Neriya98/vaccine-delivery-optimization 1c0118f4a87cf7203a0a62dcec3ebb435d967d1e \
  /data/raw/grippe-passages-aux-urgences-et-actes-sos-medecins-france.csv
# SE: Socialstyrelsen statistical database, inpatient care by main diagnosis (extract of 2024-11-25)
fetch national/se https://github.com/klpn/static-dust 90b288321438924a1c35ed1d9b852d0f8f692e6c \
  /postdata/2024-11-25-vag/data/infl_se_par_085.tsv
# PL: NIZP-PZH 'Meldunki' influenza bulletins, converted to CSV
fetch national/pl https://github.com/krzysiektr/Grypa d92abb55693bd8b46dac1e383dfb1f2f1159bb7b \
  /g_szpital_2010_2022.csv /g_szpital_2023.csv

echo "sources fetched into $W"
