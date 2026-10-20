#!/usr/bin/env bash
# Fetch the ERVISS releases that together cover every week ERVISS has published of its primary-care
# files, pinned to the ERVISS repository commit of 2026-10-09.
#
# ERVISS publishes a rolling window, and ECDC keeps each weekly release in data/snapshots/. Each window
# start was found by opening the releases on either side of every change:
#   releases up to 2024-10-18 start at 2014-W40
#   releases up to 2025-06-20 start at 2020-W25
#   releases up to 2026-05-22 start at 2021-W25
#   later releases and the current files start at 2022-W25
# The last release of each window plus the current file therefore holds the last-published value of
# every week. (The pipeline's own snapshots under data/ are the 2026-05-22 release.)
#
# Usage: bash code/09_erviss_history/fetch_releases.sh <work_dir>
set -euo pipefail
W=${1:?usage: fetch_releases.sh <work_dir>}
DIR=$W/erviss_history
STEMS="ILIARIRates sentinelTestsDetectionsPositivity nonSentinelTestsDetections"
RELEASES="2024-10-18 2025-06-20 2026-05-22"

paths=()
for s in $STEMS; do
  paths+=("/data/$s.csv")
  for d in $RELEASES; do paths+=("/data/snapshots/${d}_$s.csv"); done
done
mkdir -p "$W"
[ -d "$DIR/.git" ] || git clone -q --filter=blob:none --no-checkout https://github.com/EU-ECDC/Respiratory_viruses_weekly_data "$DIR"
git -C "$DIR" sparse-checkout set --no-cone "${paths[@]}"
git -C "$DIR" checkout -q c41550037d19ede7988bfa23ea4261abfc28b7cc
echo "ERVISS releases fetched into $DIR"
