#!/usr/bin/env bash
set -euo pipefail

echo "Do not run from Codex. Manual user execution only."

if [[ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES after reviewing official GEO files."
  exit 1
fi

# Backup PD microarray cohorts. Prefer series matrix and processed expression tables.

# GSE7621
# mkdir -p data/raw/external/gse7621_pd_sn_bulk
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7621/matrix/GSE7621_series_matrix.txt.gz" -o data/raw/external/gse7621_pd_sn_bulk/GSE7621_series_matrix.txt.gz
# shasum -a 256 data/raw/external/gse7621_pd_sn_bulk/GSE7621_series_matrix.txt.gz > data/raw/external/gse7621_pd_sn_bulk/GSE7621_series_matrix.txt.gz.sha256

# GSE8397
# mkdir -p data/raw/external/gse8397_pd_sn_bulk
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/matrix/GSE8397_series_matrix.txt.gz" -o data/raw/external/gse8397_pd_sn_bulk/GSE8397_series_matrix.txt.gz
# shasum -a 256 data/raw/external/gse8397_pd_sn_bulk/GSE8397_series_matrix.txt.gz > data/raw/external/gse8397_pd_sn_bulk/GSE8397_series_matrix.txt.gz.sha256

# GSE20292
# mkdir -p data/raw/external/gse20292_pd_sn_bulk
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20292/matrix/GSE20292_series_matrix.txt.gz" -o data/raw/external/gse20292_pd_sn_bulk/GSE20292_series_matrix.txt.gz
# shasum -a 256 data/raw/external/gse20292_pd_sn_bulk/GSE20292_series_matrix.txt.gz > data/raw/external/gse20292_pd_sn_bulk/GSE20292_series_matrix.txt.gz.sha256

echo "Review GEO filenames, uncomment only needed lightweight commands, then run manually."
