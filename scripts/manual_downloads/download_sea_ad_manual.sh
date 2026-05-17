set -euo pipefail

# DO NOT RUN FROM CODEX.
# MANUAL_HEAVY: SEA-AD processed single-nucleus RNA-seq download template.
# This script is intentionally guarded. It will not download unless you set:
#   RUN_MANUAL_DOWNLOAD=YES

echo "WARNING: This is a MANUAL_HEAVY download template for SEA-AD."
echo "WARNING: Review Allen/BKP/AWS access terms before running."
echo "WARNING: Do not run this from Codex."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Refusing to run. Set RUN_MANUAL_DOWNLOAD=YES manually to proceed."
  exit 1
fi

SEA_AD_SOURCE_URI="${SEA_AD_SOURCE_URI:-}"
if [[ -z "${SEA_AD_SOURCE_URI}" ]]; then
  echo "Set SEA_AD_SOURCE_URI to the official SEA-AD AWS/Open Data file or prefix first."
  echo "Example placeholder:"
  echo "  SEA_AD_SOURCE_URI='s3://OFFICIAL_SEA_AD_BUCKET/OFFICIAL_PROCESSED_SNRNA_FILE'"
  exit 1
fi

mkdir -p data/raw/sea_ad

echo "MANUAL_HEAVY: downloading SEA-AD from ${SEA_AD_SOURCE_URI}"
echo "Output directory: data/raw/sea_ad"

# MANUAL_HEAVY: choose one official command after filling SEA_AD_SOURCE_URI.
# aws s3 cp --no-sign-request "${SEA_AD_SOURCE_URI}" data/raw/sea_ad/
# aws s3 sync --no-sign-request "${SEA_AD_SOURCE_URI}" data/raw/sea_ad/

echo "Template complete. Uncomment exactly one MANUAL_HEAVY command after verifying the source."
