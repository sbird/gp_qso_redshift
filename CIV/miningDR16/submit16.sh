#!/bin/bash -l

# Submit directly from the HPCC terminal:
#   sbatch submit16.sh FIRST_QSO NUMBER_OF_QSOS

#SBATCH --nodes=1
#SBATCH --ntasks=1
#   SBATCH --cpus-per-task=WORKERS
#   SBATCH --mem=64gb
#SBATCH --time=72:00:00
#   SBATCH --partition=highclock
#SBATCH --output=prc_%j.out
#SBATCH --job-name=prc_dr16
#SBATCH --mail-user=rmonadi@calpoly.edu
#SBATCH --mail-type=ALL

set -euo pipefail


start=''
count=''
workers=32
cpu='highclock'
mem='64gb'
gate=1
singlet=0
injector=0;
while (($#)); do
    case "$1" in
        --start=*) start="${1#*=}" ;;
        --count=*) count="${1#*=}" ;;
        --workers=*) workers="${1#*=}" ;;
        --gate=*) gate="${1#*=}" ;;
        --singlet=*) singlet="${1#*=}" ;;
        --injector=*) injector="${1#*=}" ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done

if [[ ! "$start" =~ ^[1-9][0-9]*$ || ! "$count" =~ ^[1-9][0-9]*$ ]]; then
    echo 'Required options: --start=POSITIVE_INTEGER --count=POSITIVE_INTEGER' >&2
    echo 'Example: sbatch submit16.sh --start=1 --count=5000 --workers=32 --gate=1 --singlet=0' >&2
    exit 2
fi

export CIV_DOUBLET_COVERAGE_GATE="$gate"
export CIV_SINGLET_MODEL="$singlet"
export INJECTOR="$injector"

if [[ "$gate" == 1 ]]; then gate_text=ON; else gate_text=OFF; fi
if [[ "$singlet" == 1 ]]; then singlet_text=ON; else singlet_text=OFF; fi
echo "Doublet coverage gate: $gate_text"
echo "Singlet model competition: $singlet_text"
echo "Injector: $injector"
module load matlab
matlab -batch "ind_S=$start; num_quasars=$count; cores = $workers;  injector = $injector; run_script_dr16_fixed; exit"
