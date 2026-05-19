#!/usr/bin/env bash
set -euo pipefail

echo "Do not run from Codex. Manual user execution only."

if [[ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES after reviewing official GEO files."
  exit 1
fi

RAW_DIR="data/raw/external/gse20186_pd_superseries"
mkdir -p "${RAW_DIR}"

# GEO accession: GSE20186
# Prefer series/subseries matrix files, SOFT/MINiML metadata, platform annotation, and processed expression tables.
# Series matrix example:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20186/matrix/GSE20186_series_matrix.txt.gz" -o "${RAW_DIR}/GSE20186_series_matrix.txt.gz"
# wget -O "${RAW_DIR}/GSE20186_series_matrix.txt.gz" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20186/matrix/GSE20186_series_matrix.txt.gz"
# SOFT:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20186/soft/GSE20186_family.soft.gz" -o "${RAW_DIR}/GSE20186_family.soft.gz"
# MINiML:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20186/miniml/GSE20186_family.xml.tgz" -o "${RAW_DIR}/GSE20186_family.xml.tgz"
# Supplementary directory:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20186/suppl/" -o "${RAW_DIR}/GSE20186_suppl_listing.html"
# Checksum placeholder:
# shasum -a 256 "${RAW_DIR}/GSE20186_series_matrix.txt.gz" > "${RAW_DIR}/GSE20186_series_matrix.txt.gz.sha256"

echo "Review GEO filenames, uncomment only needed lightweight commands, then run manually."
