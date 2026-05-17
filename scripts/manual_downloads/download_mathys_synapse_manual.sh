set -euo pipefail

# DO NOT RUN FROM CODEX.
# MANUAL_HEAVY: Mathys 2019 AD single-cell/single-nucleus download template.
# This script is intentionally guarded. It will not download unless you set:
#   RUN_MANUAL_DOWNLOAD=YES

echo "WARNING: This is a MANUAL_HEAVY download template for Mathys 2019 / GSE138852 / Synapse-linked data."
echo "WARNING: Confirm access permissions, terms, and correct file IDs before running."
echo "WARNING: Do not run this from Codex."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Refusing to run. Set RUN_MANUAL_DOWNLOAD=YES manually to proceed."
  exit 1
fi

MATHYS_SYNAPSE_ID="${MATHYS_SYNAPSE_ID:-}"
if [[ -z "${MATHYS_SYNAPSE_ID}" ]]; then
  echo "Set MATHYS_SYNAPSE_ID to the approved Synapse file/project ID first."
  echo "GEO accession context: GSE138852."
  echo "Example placeholder:"
  echo "  MATHYS_SYNAPSE_ID='synXXXXXXXX'"
  exit 1
fi

mkdir -p data/raw/mathys_2019

echo "MANUAL_HEAVY: downloading Mathys 2019 resource ${MATHYS_SYNAPSE_ID}"
echo "Output directory: data/raw/mathys_2019"

# MANUAL_HEAVY: uncomment after Synapse login/access approval and ID verification.
# synapse get "${MATHYS_SYNAPSE_ID}" --downloadLocation data/raw/mathys_2019

echo "Template complete. Uncomment the MANUAL_HEAVY command after verifying access."
