#!/usr/bin/env bash
# Rebuild the ERVISS primary-care series merged across releases, from 2014/15:
#   data/primary_care_ili_rates/       ILI and ARI rates, coverage, comparison with RespiCompass ILI+, figure
#   data/primary_care_flu_positivity/  influenza tests and detections (sentinel, non-sentinel), coverage, figure
# Run from the repository root:  bash code/09_erviss_history/build_all.sh <work_dir>
# <work_dir> receives about 70 MB of ERVISS releases. Needs git, Python 3 (standard library only) and R
# from the project's renv (ggplot2, dplyr, tidyr, ISOweek).
set -euo pipefail
W=${1:?usage: build_all.sh <work_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)

bash "$HERE/fetch_releases.sh" "$W"
python3 -I "$HERE/merge_releases.py" "$W/erviss_history/data" data
python3 -I "$HERE/compare_iliplus.py" data
Rscript "$HERE/plot_history.R" data
