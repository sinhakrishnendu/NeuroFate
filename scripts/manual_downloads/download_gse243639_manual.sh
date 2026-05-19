#!/usr/bin/env bash
set -euo pipefail

echo "Do not run from Codex. Manual user execution only."
echo "Dataset: GSE243639 PD substantia nigra pars compacta snRNA-seq."

if [[ "${RUN_MANUAL_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set RUN_MANUAL_DOWNLOAD=YES after editing official filenames and checksums."
  exit 1
fi

OUTDIR="data/raw/external/gse243639_pd_snpc"
mkdir -p "${OUTDIR}"

echo "Expected local path: ${OUTDIR}"
echo "Edit this script with official GEO supplementary filenames before running."

# MANUAL_HEAVY checksum placeholders:
# EXPECTED_SHA256_COUNTS="PASTE_OFFICIAL_OR_USER_COMPUTED_SHA256"
# EXPECTED_SHA256_METADATA="PASTE_OFFICIAL_OR_USER_COMPUTED_SHA256"

# MANUAL_HEAVY examples only; keep commented until official filenames are verified:
# wget -O "${OUTDIR}/OFFICIAL_COUNTS_FILE" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243639/suppl/OFFICIAL_COUNTS_FILE"
# wget -O "${OUTDIR}/OFFICIAL_METADATA_FILE" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243639/suppl/OFFICIAL_METADATA_FILE"
# shasum -a 256 "${OUTDIR}/OFFICIAL_COUNTS_FILE"
# shasum -a 256 "${OUTDIR}/OFFICIAL_METADATA_FILE"
