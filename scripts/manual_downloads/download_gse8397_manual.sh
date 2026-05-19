#!/usr/bin/env bash
set -euo pipefail

echo "Do not run from Codex. Manual user execution only."

if [[ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES after reviewing official GEO files."
  exit 1
fi

RAW_DIR="data/raw/external/gse8397_pd_brain_regions"
mkdir -p "${RAW_DIR}"

# GEO accession: GSE8397
# Prefer series matrix, SOFT/MINiML, platform annotation, and processed expression tables.
# Series matrix:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/matrix/GSE8397_series_matrix.txt.gz" -o "${RAW_DIR}/GSE8397_series_matrix.txt.gz"
# SOFT:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/soft/GSE8397_family.soft.gz" -o "${RAW_DIR}/GSE8397_family.soft.gz"
# MINiML:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/miniml/GSE8397_family.xml.tgz" -o "${RAW_DIR}/GSE8397_family.xml.tgz"
# Checksum placeholder:
# shasum -a 256 "${RAW_DIR}/GSE8397_series_matrix.txt.gz" > "${RAW_DIR}/GSE8397_series_matrix.txt.gz.sha256"

echo "Review official filenames, uncomment only lightweight metadata/processed commands, then run manually."
