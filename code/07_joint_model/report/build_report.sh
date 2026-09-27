#!/bin/sh
# build_report.sh -- the printable A4 report: REPORT.md -> HTML (pandoc) -> PDF (headless Chromium).
# Run report_figures.R first. Writes documentation/joint_model_report.pdf. The fonts are Source Serif 4
# and Source Sans 3 (Google Fonts); without them the Liberation faces stand in.
# CHROME may point at any Chrome or Chromium binary.
set -e
cd "$(dirname "$0")/../../.."
OUT=output/joint_model/report
CHROME=${CHROME:-$(ls /opt/pw-browsers/chromium-*/chrome-linux/chrome 2>/dev/null | head -1)}
cp code/07_joint_model/report/report.css "$OUT/"
pandoc code/07_joint_model/report/REPORT.md -f markdown -t html5 \
  --template code/07_joint_model/report/template.html --css report.css \
  --metadata pagetitle="Joint EU/EEA influenza model: status report" -o "$OUT/report.html"
"$CHROME" --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=documentation/joint_model_report.pdf "file://$PWD/$OUT/report.html" 2>/dev/null
echo "wrote documentation/joint_model_report.pdf"
