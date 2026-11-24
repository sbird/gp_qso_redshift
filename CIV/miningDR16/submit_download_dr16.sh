#!/usr/bin/env bash
#SBATCH --job-name=download-dr16
#SBATCH --output=download-dr16-%j.out
#SBATCH --error=download-dr16-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1

# Submit from the repository root:
#   sbatch submit_download_dr16.sh
#
# Optional overrides:
#   sbatch --export=ALL,LIST_FILE=/path/to/dr16Spec.list,DESTINATION=/path/to/data/dr16/sas/dr16 submit_download_dr16.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIST_FILE="${LIST_FILE:-${REPO_DIR}/data/dr16Spec.list}"
DESTINATION="${DESTINATION:-${REPO_DIR}/data/dr16/sas/dr16}"

if [[ ! -f "${LIST_FILE}" ]]; then
    echo "List file not found: ${LIST_FILE}" >&2
    exit 1
fi

mkdir -p "${DESTINATION}"
cd "${DESTINATION}"

# The list contains full https://data.sdss.org/sas/dr16/... URLs.  Removing
# the first two remote directories (sas/dr16) makes files land at
# DESTINATION/eboss/... and DESTINATION/sdss/..., matching the MATLAB loaders.
wget --continue --force-directories --no-host-directories --cut-dirs=2 \
    --timeout=60 --tries=5 --waitretry=10 --input-file="${LIST_FILE}"
