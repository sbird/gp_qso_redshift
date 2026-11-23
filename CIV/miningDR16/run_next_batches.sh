#!/usr/bin/env bash
set -euo pipefail

# First new batch waits for your existing final job.
dependency=(--dependency=afterok:31563414 --kill-on-invalid-dep=yes)

for start in 30001 35001 40001 45001 50001; do
    echo "Submitting batch $start–$((start + 4999))"

	if sbatch --wait "${dependency[@]}" submit16.sh "$start" 5000; then
		echo "Batch starting at $start completed successfully."
	else
		echo "Batch starting at $start failed; stopping." >&2
		exit 1
	fi

# Subsequent submissions happen only after sbatch --wait returns.
	dependency=()
done

