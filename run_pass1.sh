#!/usr/bin/env bash
###############################################################################
# run_pass1.sh -- QC, integration/clustering, BoneMarrowMap (scripts 01-03)
#
# Start it so it survives a dropped SSH session:
#   cd /home/danielpy/data/scripts
#   nohup bash run_pass1.sh > ../Output/logs/run_pass1.out 2>&1 &
#   tail -f ../Output/logs/run_pass1.out          # watch progress
#
# Each script writes its own log to Output/logs/<script>.log.
# If a script fails, the pipeline stops there; fix the problem and rerun.
# Already-finished steps can be skipped by commenting them out below.
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${FPD_BASE_DIR:-/home/danielpy/data}"
LOG_DIR="${BASE_DIR}/Output/logs"
mkdir -p "${LOG_DIR}"

# Use GNU time (if installed) to record each script's peak memory
TIMER=""
if [ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1; then
  TIMER="/usr/bin/time -v"
fi

run_step () {
  local script="$1"
  local log="${LOG_DIR}/${script%.R}.log"
  echo "[$(date '+%F %T')] START ${script}  (log: ${log})"
  if ${TIMER} Rscript "${SCRIPT_DIR}/${script}" > "${log}" 2>&1; then
    echo "[$(date '+%F %T')] DONE  ${script}"
    grep -E "Maximum resident set size|Elapsed \(wall clock\)" "${log}" || true
  else
    echo "[$(date '+%F %T')] FAILED ${script} -- last lines of ${log}:"
    tail -n 25 "${log}"
    exit 1
  fi
}

run_step 01_sample_sheet_qc.R
run_step 02_integrate_cluster.R
run_step 03_bonemarrowmap.R

echo "[$(date '+%F %T')] Pass 1 finished."
echo "Next: review Output/, then copy metadata/downsample_overrides_TEMPLATE.csv"
echo "      to metadata/downsample_overrides.csv, edit CellsPerDonor, and run run_pass2.sh"

# ---- Optional: upload results to S3 (fill in your bucket, then uncomment) ----
# aws s3 sync "${BASE_DIR}/Output" "s3://YOUR-BUCKET/YOUR-PREFIX/Output" --exclude "checkpoints/*"
