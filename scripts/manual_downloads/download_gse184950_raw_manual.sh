#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: GSE184950 RAW archive download template."
echo "Do not run from Codex. Manual user execution only."
echo "This script exits unless RUN_MANUAL_DOWNLOAD=YES."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES only after reviewing GEO GSE184950."
  exit 1
fi

RAW_DIR="data/raw/external/gse184950_pd_sn"
mkdir -p "${RAW_DIR}"

echo "Review GEO and checksum expectations before uncommenting commands."
# GEO page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184950
# Expected local file:
# data/raw/external/gse184950_pd_sn/GSE184950_RAW.tar
#
# MANUAL_HEAVY download templates. Keep commented until manually reviewed:
# curl -L -o "${RAW_DIR}/GSE184950_RAW.tar" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184950/suppl/GSE184950_RAW.tar"
# wget -O "${RAW_DIR}/GSE184950_RAW.tar" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184950/suppl/GSE184950_RAW.tar"
# shasum -a 256 "${RAW_DIR}/GSE184950_RAW.tar" > metadata/checksums/gse184950_raw_sha256.tsv

echo "Template complete. No download command is active by default."
