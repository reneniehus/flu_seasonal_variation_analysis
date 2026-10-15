#!/usr/bin/env bash
# Rebuild data/flu_hospitalisation_indicators/ from the pinned public sources.
# Run from the repository root:  bash code/08_flu_hospitalisation/build_all.sh <work_dir>
# <work_dir> receives the downloaded sources (about 150 MB); it can live outside the repository.
# Needs git, Python 3 (standard library only) and R with the renv library (ggplot2, dplyr, ISOweek).
set -euo pipefail
W=${1:?usage: build_all.sh <work_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=data/flu_hospitalisation_indicators

bash "$HERE/fetch_sources.sh" "$W"
python3 -I "$HERE/build_erviss.py" "$W/erviss/data" "$OUT"
python3 -I "$HERE/build_national.py" "$W" "$OUT"
python3 -I "$HERE/build_coverage.py" "$OUT"
python3 -I "$HERE/reconcile.py" "$OUT"
Rscript "$HERE/plot_coverage.R" "$OUT"
Rscript "$HERE/plot_series.R" "$OUT"
