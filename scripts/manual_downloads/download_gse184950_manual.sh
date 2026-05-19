#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: GSE184950 acquisition template."
echo "Do not run from Codex. Manual user execution only."
echo "This script is guarded and exits unless RUN_MANUAL_DOWNLOAD=YES."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES only when you are ready to manually acquire GSE184950."
  exit 1
fi

RAW_DIR="data/raw/external/gse184950_pd_sn"
mkdir -p "${RAW_DIR}"

echo "Review GEO GSE184950 manually before uncommenting any command."
# GEO page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184950
# Candidate supplementary files may include:
# - GSE184950_RAW.tar
# - GSE184950_add2.xlsx
#
# MANUAL_HEAVY examples, keep commented until official filenames are verified:
# wget -P "${RAW_DIR}" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184950/suppl/GSE184950_RAW.tar"
# wget -P "${RAW_DIR}" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184950/suppl/GSE184950_add2.xlsx"
# shasum -a 256 "${RAW_DIR}"/* > metadata/checksums/gse184950_sha256.tsv

echo "Template complete. No command is uncommented by default."
