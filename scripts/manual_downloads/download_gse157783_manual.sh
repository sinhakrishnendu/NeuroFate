#!/usr/bin/env bash
set -euo pipefail

echo "Do not run from Codex. Manual user execution only."

if [[ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES after reviewing official GEO files."
  exit 1
fi

RAW_DIR="data/raw/external/gse157783_pd_midbrain_snrna_optional"
mkdir -p "${RAW_DIR}"

# GEO accession: GSE157783
# Optional small snRNA cohort. Use only processed matrices if available.
# Do not run SRA tools or FASTQ processing for NeuroFate Phase 33.
# Series matrix example:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/matrix/GSE157783_series_matrix.txt.gz" -o "${RAW_DIR}/GSE157783_series_matrix.txt.gz"
# Supplementary directory:
# curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/suppl/" -o "${RAW_DIR}/GSE157783_suppl_listing.html"
# Checksum placeholder:
# shasum -a 256 "${RAW_DIR}/GSE157783_series_matrix.txt.gz" > "${RAW_DIR}/GSE157783_series_matrix.txt.gz.sha256"

echo "Review GEO filenames, uncomment only needed lightweight commands, then run manually."
