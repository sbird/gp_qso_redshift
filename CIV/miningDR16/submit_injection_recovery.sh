#!/bin/bash -l
# Submit from the project root, for example:
#   sbatch submit_injection_recovery.sh 1 25000 5000
# Arguments: first filtered-test row, last filtered-test row, injection count,
# series number, minimum absorber redshift, minimum injected REW1548, and minimum
# parent-spectrum S/N. Series 2 excludes saved series 1 sightlines.
# Override resources with sbatch options if your cluster requires different limits.
#SBATCH --job-name=civ_inject
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --partition=epyc
#SBATCH --output=civ_inject_%j.out

set -euo pipefail
if [[ $# -lt 2 || $# -gt 8 ]]; then
    echo "Usage: sbatch $0 FIRST_TEST_ROW LAST_TEST_ROW [INJECTION_COUNT] [SERIES] [MIN_ZABS] [MIN_REW1548] [MIN_SNR] [RESUME]" >&2
    exit 2
fi
for value in "${@:1:4}"; do
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo 'All arguments must be positive integers without leading zeros.' >&2
        exit 2
    fi
done
if [[ $# -ge 5 && ! "$5" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo 'MIN_ZABS must be a nonnegative number.' >&2
    exit 2
fi
if [[ $# -ge 6 && ! "$6" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo 'MIN_REW1548 must be a nonnegative number.' >&2
    exit 2
fi
if [[ $# -ge 7 && ! "$7" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo 'MIN_SNR must be a nonnegative number.' >&2
    exit 2
fi
if [[ $# -ge 8 && "$8" != 1 ]]; then
    echo 'RESUME, when supplied, must be 1.' >&2
    exit 2
fi
if (( $1 > $2 )); then
    echo 'FIRST_TEST_ROW must not exceed LAST_TEST_ROW.' >&2
    exit 2
fi
export CIV_INJECTION_START="$1"
export CIV_INJECTION_END="$2"
export CIV_INJECTION_COUNT="${3:-5000}"
export CIV_INJECTION_SERIES="${4:-1}"
export CIV_INJECTION_MIN_ZABS="${5:-}"
export CIV_INJECTION_MIN_REW1548="${6:-}"
export CIV_INJECTION_SNR_THRESHOLD="${7:-}"
export CIV_INJECTION_RESUME="${8:-}"
coverage_gate="${CIV_DOUBLET_COVERAGE_GATE:-1}"
singlet_model="${CIV_SINGLET_MODEL:-0}"
if [[ "$coverage_gate" != 0 && "$coverage_gate" != 1 ]]; then
    echo 'CIV_DOUBLET_COVERAGE_GATE must be 0 or 1.' >&2
    exit 2
fi
if [[ "$singlet_model" != 0 && "$singlet_model" != 1 ]]; then
    echo 'CIV_SINGLET_MODEL must be 0 or 1.' >&2
    exit 2
fi
export CIV_DOUBLET_COVERAGE_GATE="$coverage_gate"
export CIV_SINGLET_MODEL="$singlet_model"
echo "Doublet coverage gate: $coverage_gate; singlet model competition: $singlet_model."
module load matlab
matlab -batch "run_script_dr16_fixed"
