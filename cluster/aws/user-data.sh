#!/bin/bash
## =====================================================================
##  The validation on one AWS EC2 instance (Amazon Linux 2023), as user data
##
##  Runs the validation in the published image, copies the log and the
##  results to S3, and shuts the instance down; launched with shutdown
##  behaviour "terminate", it then deletes itself. A safety cutoff shuts
##  it down after MAX_HOURS whatever happens. cluster/README.md, "On AWS".
##
##  Set per run (the launch command fills them in):
##    BUCKET      the S3 bucket for the log and the results
##    RUN_NAME    a folder in the bucket for this run
##    SMOKE       1 for the smoke test, 0 for a real run
##    STUDIES     the studies to run ("" for all), as for run_all.R
##    MAX_HOURS   the safety cutoff
## =====================================================================
BUCKET="__BUCKET__"
RUN_NAME="__RUN_NAME__"
SMOKE="__SMOKE__"
STUDIES="__STUDIES__"
MAX_HOURS="__MAX_HOURS__"

IMAGE="ghcr.io/dm-k1/flatcor-validation:main"
LOG=/var/log/flatcor-run.log
exec > "$LOG" 2>&1
set -x

## the safety cutoff, first of all
shutdown -h "+$((MAX_HOURS * 60))"

## whatever happens below: copy the log and the results to S3, then shut down
finish() {
  cd /opt/flatcor-validation 2>/dev/null || true
  results=results
  [ "$SMOKE" = "1" ] && results=cache/smoke_results
  aws s3 sync "$results" "s3://$BUCKET/$RUN_NAME/results" --only-show-errors || true
  aws s3 sync cache/logs "s3://$BUCKET/$RUN_NAME/logs" --only-show-errors || true
  aws s3 cp "$LOG" "s3://$BUCKET/$RUN_NAME/run.log" --only-show-errors || true
  shutdown -h now
}
trap finish EXIT

## Docker, git and the project
dnf install -y docker git
systemctl start docker
git clone https://github.com/dm-k1/flatcor-validation.git /opt/flatcor-validation
cd /opt/flatcor-validation || exit 1
git log --oneline -1
docker pull "$IMAGE"

## the run
echo "cores: $(nproc)"
# shellcheck disable=SC2086  # STUDIES is a list of words on purpose
docker run --rm -v "$PWD":/project -e VALIDATION_SMOKE="$SMOKE" "$IMAGE" Rscript run_all.R $STUDIES
echo "run_all.R exited with status $?"
