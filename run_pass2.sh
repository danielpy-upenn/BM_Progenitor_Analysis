#!/usr/bin/env bash
###############################################################################
# run_pass2.sh -- differential expression (script 04)
#
# Requires metadata/downsample_overrides.csv (copied from the TEMPLATE that
# script 03 wrote, with CellsPerDonor edited).
#
#   cd /home/danielpy/data/scripts
#   nohup bash run_pass2.sh > ../Output/logs/run_pass2.out 2>&1 &
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${FPD_BASE_DIR:-/home/danielpy/data}"
LOG_DIR="${BASE_DIR}/Output/logs"
mkdir -p "${LOG_DIR}"

TIMER=""
if [ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1; then
  TIMER="/usr/bin/time -v"
fi

script=04_differential_expression.R
log="${LOG_DIR}/${script%.R}.log"
echo "[$(date '+%F %T')] START ${script}  (log: ${log})"
if ${TIMER} Rscript "${SCRIPT_DIR}/${script}" > "${log}" 2>&1; then
  echo "[$(date '+%F %T')] DONE  ${script}"
  grep -E "Maximum resident set size|Elapsed \(wall clock\)" "${log}" || true
else
  echo "[$(date '+%F %T')] FAILED ${script} -- last lines of ${log}:"
  tail -n 25 "${log}"
  exit 1
fi

# ---- Optional: upload results to S3 (fill in your bucket, then uncomment) ----
# aws s3 sync "${BASE_DIR}/Output" "s3://YOUR-BUCKET/YOUR-PREFIX/Output" --exclude "checkpoints/*"
