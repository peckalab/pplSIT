#!/usr/bin/env bash
# Dry-run by default. Uses the existing pipe-dev environment.
set -euo pipefail
batch_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$batch_dir/../../../.." && pwd)"
stage="${1:-}"
mode="${2:---dry-run}"
case "$stage" in
  pilot-active) configs=("$batch_dir/active_sort_qc.yaml" "$batch_dir/pilot_active.yaml") ;;
  pilot-passive) configs=("$batch_dir/passive_reuse_qc.yaml" "$batch_dir/pilot_passive.yaml") ;;
  active) configs=("$batch_dir/active_sort_qc.yaml") ;;
  active-sync) configs=("$batch_dir/active_sort_qc.yaml" "$batch_dir/active_sync.yaml") ;;
  passive) configs=("$batch_dir/passive_reuse_qc.yaml") ;;
  *) echo 'Usage: bash run.sh {pilot-active|pilot-passive|active|active-sync|passive} [--dry-run|--execute]' >&2; exit 2 ;;
esac
case "$mode" in
  --dry-run) options=(--dry-run --nolock) ;;
  --execute) options=(--rerun-incomplete) ;;
  *) echo 'Expected --dry-run or --execute' >&2; exit 2 ;;
esac
if (( $# > 2 )); then echo 'Unexpected extra arguments' >&2; exit 2; fi
export PATH="/opt/miniconda3/envs/pipe-dev/bin:/opt/miniconda3/bin:$PATH"
if [[ "$mode" == --execute && "$stage" == *active ]]; then
  nvidia-smi --id=3,1 --query-gpu=index,memory.free,utilization.gpu --format=csv
fi
cd "$repo_dir/workflow"
exec snakemake all --snakefile Snakefile --configfile "${configs[@]}" \
  --use-conda --cores 64 \
  --set-resources do_kilosort_stream:_cores=16 \
  --resources kilosort_gpu=2 bombcell=1 io_heavy=1 \
  "${options[@]}"
