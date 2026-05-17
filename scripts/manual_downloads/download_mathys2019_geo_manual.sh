#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: Mathys 2019 / GSE138852 acquisition template."
echo "DO NOT RUN FROM CODEX. Run manually only after reviewing GEO/Synapse access terms."
echo "No command in this file should download data unless RUN_MANUAL_DOWNLOAD=YES is set."

: "${RUN_MANUAL_DOWNLOAD:=NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "RUN_MANUAL_DOWNLOAD is not YES. Exiting without download."
  echo "Review and copy the commented wget/curl templates below if you want to download manually."
  exit 0
fi

mkdir -p data/raw/external/mathys_2019
mkdir -p metadata/checksums

echo "GEO landing page:"
echo "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852"
echo "GEO supplementary FTP directory:"
echo "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/"

echo "MANUAL_HEAVY templates are intentionally commented:"

# MANUAL_HEAVY - download GEO supplementary archive after reviewing file names and terms.
# curl -L \
#   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/" \
#   -o data/raw/external/mathys_2019/GSE138852_supplementary_listing.html

# MANUAL_HEAVY - example only; replace FILE_NAME with the actual GEO supplementary file.
# wget -O data/raw/external/mathys_2019/FILE_NAME \
#   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/FILE_NAME"

# MANUAL_HEAVY - example only; replace FILE_NAME with the actual GEO supplementary file.
# curl -L \
#   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/FILE_NAME" \
#   -o data/raw/external/mathys_2019/FILE_NAME

# LIGHTWEIGHT CHECKSUM PLACEHOLDERS - run manually after download.
# md5 data/raw/external/mathys_2019/FILE_NAME > metadata/checksums/mathys2019_FILE_NAME.md5
# shasum -a 256 data/raw/external/mathys_2019/FILE_NAME > metadata/checksums/mathys2019_FILE_NAME.sha256

echo "Expected local paths:"
echo "  data/raw/external/mathys_2019/"
echo "  data/interim/external/mathys_2019/"
echo "  results/logs/36_inspect_mathys2019.log"
echo "  results/tables/mathys2019_metadata_summary.tsv"
echo "  data/interim/external/mathys_2019/mathys_var_genes.tsv"
