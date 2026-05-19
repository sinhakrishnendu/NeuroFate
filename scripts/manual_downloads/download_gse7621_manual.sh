#!/usr/bin/env bash
set -euo pipefail

echo "Do not run from Codex. Manual user execution only."

if [[ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES after reviewing official GEO files."
  exit 1
fi

RAW_DIR="data/raw/external/gse7621_pd_sn_bulk"
mkdir -p "${RAW_DIR}"

# GEO accession: GSE7621
# Prefer series matrix, SOFT/MINiML, platform annotation, and processed expression tables.
# Series matrix:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7621/matrix/GSE7621_series_matrix.txt.gz" -o "${RAW_DIR}/GSE7621_series_matrix.txt.gz"
# Platform annotation, only if the parsed platform is not already mapped locally:
# PLATFORM_ID="<GPL_FROM_PHASE37_PLATFORM_SUMMARY>"
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/platforms/${PLATFORM_ID:0:5}nnn/${PLATFORM_ID}/annot/${PLATFORM_ID}.annot.gz" -o "data/raw/platforms/${PLATFORM_ID}/${PLATFORM_ID}.annot.gz"
# SOFT:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7621/soft/GSE7621_family.soft.gz" -o "${RAW_DIR}/GSE7621_family.soft.gz"
# MINiML:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7621/miniml/GSE7621_family.xml.tgz" -o "${RAW_DIR}/GSE7621_family.xml.tgz"
# Checksum placeholder:
# shasum -a 256 "${RAW_DIR}/GSE7621_series_matrix.txt.gz" > "${RAW_DIR}/GSE7621_series_matrix.txt.gz.sha256"
# shasum -a 256 "data/raw/platforms/${PLATFORM_ID}/${PLATFORM_ID}.annot.gz" > "data/raw/platforms/${PLATFORM_ID}/${PLATFORM_ID}.annot.gz.sha256"

echo "Review official filenames, uncomment only lightweight metadata/processed commands, then run manually."
