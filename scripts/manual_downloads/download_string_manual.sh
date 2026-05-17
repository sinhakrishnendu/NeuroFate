set -euo pipefail

# DO NOT RUN FROM CODEX.
# MANUAL_HEAVY: STRING human protein interaction links download template.
# This script is intentionally guarded. It will not download unless you set:
#   RUN_MANUAL_DOWNLOAD=YES

echo "WARNING: This is a MANUAL_HEAVY download template for STRING."
echo "WARNING: Review STRING license and citation requirements before running."
echo "WARNING: Do not run this from Codex."

RUN_MANUAL_DOWNLOAD="${RUN_MANUAL_DOWNLOAD:-NO}"
if [[ "${RUN_MANUAL_DOWNLOAD}" != "YES" ]]; then
  echo "Refusing to run. Set RUN_MANUAL_DOWNLOAD=YES manually to proceed."
  exit 1
fi

STRING_SOURCE_URL="${STRING_SOURCE_URL:-}"
if [[ -z "${STRING_SOURCE_URL}" ]]; then
  echo "Set STRING_SOURCE_URL to the official STRING human protein links file URL first."
  echo "Example placeholder:"
  echo "  STRING_SOURCE_URL='https://string-db.org/path/to/official/human.protein.links.txt.gz'"
  exit 1
fi

mkdir -p data/external/string

echo "MANUAL_HEAVY: downloading STRING human interaction links from ${STRING_SOURCE_URL}"
echo "Output directory: data/external/string"

# MANUAL_HEAVY: uncomment after verifying the official STRING URL.
# curl -L "${STRING_SOURCE_URL}" -o data/external/string/string_human_protein_links.txt.gz

echo "Template complete. Uncomment the MANUAL_HEAVY command after verifying the source."
