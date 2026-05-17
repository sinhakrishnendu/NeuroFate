set -euo pipefail

# DO NOT RUN FROM CODEX.
# MANUAL_HEAVY: ROSMAP / AD Knowledge Portal / Synapse download template.
# This script is intentionally guarded. It will not download unless you set:
#   RUN_MANUAL_DOWNLOAD=YES

echo "WARNING: This is a MANUAL_HEAVY download template for ROSMAP controlled-access data."
echo "WARNING: Complete DUC/Synapse permissions and verify file IDs before running."
echo "WARNING: Do not run this from Codex."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Refusing to run. Set RUN_MANUAL_DOWNLOAD=YES manually to proceed."
  exit 1
fi

ROSMAP_SYNAPSE_ID="${ROSMAP_SYNAPSE_ID:-}"
if [[ -z "${ROSMAP_SYNAPSE_ID}" ]]; then
  echo "Set ROSMAP_SYNAPSE_ID to the approved Synapse file/project ID first."
  echo "Example placeholder:"
  echo "  ROSMAP_SYNAPSE_ID='synXXXXXXXX'"
  exit 1
fi

mkdir -p data/raw/rosmap

echo "MANUAL_HEAVY: downloading ROSMAP resource ${ROSMAP_SYNAPSE_ID}"
echo "Output directory: data/raw/rosmap"

# MANUAL_HEAVY: uncomment after DUC/Synapse approval and ID verification.
# synapse get "${ROSMAP_SYNAPSE_ID}" --downloadLocation data/raw/rosmap

echo "Template complete. Uncomment the MANUAL_HEAVY command after verifying access."
