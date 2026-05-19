#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: GSE174367 AD replication acquisition template."
echo "Do not run from Codex. Manual user execution only."
echo "All download commands remain commented until the user reviews GEO files."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES only after reviewing GEO GSE174367 and data-use terms."
  exit 1
fi

RAW_DIR="data/raw/external/gse174367_ad_multiomics"
mkdir -p "${RAW_DIR}"

# GEO page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367
# GEO supplementary FTP directory:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/
#
# MANUAL_HEAVY examples only. Keep commented until official filenames are selected:
# curl -L -o "${RAW_DIR}/GSE174367_series_matrix.txt.gz" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/matrix/GSE174367_series_matrix.txt.gz"
# wget -P "${RAW_DIR}" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/<reviewed_file_name>"
# shasum -a 256 "${RAW_DIR}"/* > metadata/checksums/gse174367_sha256.tsv

echo "Template complete. No active download command is included."
