#!/bin/bash
## =====================================================================
##  The validation on one AWS EC2 instance (Amazon Linux 2023), as user data
##
##  Runs the validation in the published image, copies the log and the
##  results to S3, and shuts the instance down; launched with shutdown
##  behaviour "terminate", it then deletes itself. A safety cutoff shuts
##  it down after MAX_HOURS whatever happens. cluster/README.md, "On AWS".
##
##  A spot instance can be taken back by AWS at any moment, and its disk
##  with it. So a real run copies the pipeline's store (every finished
##  piece) to S3 every 15 minutes, and starts by copying back what an
##  earlier instance with the same RUN_NAME left there: launched again
##  after an interruption, it continues where that one was, losing at most
##  15 minutes of work. (The smoke test always starts from scratch.)
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
  [ "$SMOKE" = "1" ] || backup_store
  aws s3 sync "$results" "s3://$BUCKET/$RUN_NAME/results" --only-show-errors || true
  aws s3 sync cache/logs "s3://$BUCKET/$RUN_NAME/logs" --only-show-errors || true
  aws s3 cp "$LOG" "s3://$BUCKET/$RUN_NAME/run.log" --only-show-errors || true
  shutdown -h now
}
trap finish EXIT

## the store to S3 (targets writes each piece to a temporary file first and
## then moves it in place, so a copy never holds half a piece)
backup_store() {
  aws s3 sync _targets "s3://$BUCKET/$RUN_NAME/store" \
    --exclude "scratch/*" --only-show-errors || true
}

## Docker, git and the project
dnf install -y docker git
systemctl start docker
git clone https://github.com/dm-k1/flatcor-validation.git /opt/flatcor-validation
cd /opt/flatcor-validation || exit 1
git log --oneline -1
docker pull "$IMAGE"

## a real run: first what an earlier instance of this run left, then a copy
## every 15 minutes while it runs
if [ "$SMOKE" != "1" ]; then
  aws s3 sync "s3://$BUCKET/$RUN_NAME/store" _targets --only-show-errors
  ls _targets/meta 2>/dev/null && echo "continuing from the store in S3"
  while true; do sleep 900; backup_store; done &
fi

## the run
echo "cores: $(nproc)"
# shellcheck disable=SC2086  # STUDIES is a list of words on purpose
docker run --rm -v "$PWD":/project -e VALIDATION_SMOKE="$SMOKE" "$IMAGE" Rscript run_all.R $STUDIES
echo "run_all.R exited with status $?"
