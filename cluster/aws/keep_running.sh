#!/bin/bash
## =====================================================================
##  Keeps one AWS run going while AWS takes spot machines back
##
##    keep_running.sh <bucket> <run name> <instance id> <fleet file>
##
##  Run on your own computer, in the background, while a run is on AWS
##  (it needs the computer awake and the AWS login valid). Every 2 minutes
##  it looks at the run's machine. While it runs, nothing happens. When it
##  is gone, the end of the run's log in S3 (s3://<bucket>/<run name>/run.log,
##  written by user-data.sh) says why:
##
##    "run_all.R exited with status 0"      the run finished: stop
##    status 1 and R's "Execution halted"   run_all.R stopped on an error in
##                                          the code: stop and show the log,
##                                          since a new machine would only
##                                          stop on the same error
##    anything else, or no status           AWS took the machine back (or it
##                                          stopped another way): launch a
##                                          new one from <fleet file> (an EC2
##                                          Fleet request with the same user
##                                          data, so the same RUN_NAME: it
##                                          continues from the store in S3)
##
##  At most 12 new machines. Regions: BUCKET_REGION (the bucket's) and
##  MACHINE_REGION (where the fleet launches), both by default
##  AWS_DEFAULT_REGION; AWS_PROFILE as usual. checks/AWS_RUNBOOK.md (local)
##  has the rest.
## =====================================================================
bucket=$1; run=$2; id=$3; fleet=$4
if [ -z "$bucket" ] || [ -z "$run" ] || [ -z "$id" ] || [ ! -f "$fleet" ]; then
  echo "usage: keep_running.sh <bucket> <run name> <instance id> <fleet file>" >&2
  exit 2
fi
bucket_region=${BUCKET_REGION:-$AWS_DEFAULT_REGION}
machine_region=${MACHINE_REGION:-$AWS_DEFAULT_REGION}
wait_seconds=${KEEPER_WAIT_SECONDS:-120}  # shorter only when testing
relaunches=0

run_log() {
  aws s3 cp --region "$bucket_region" "s3://$bucket/$run/run.log" - 2>/dev/null
}

## the status run_all.R ended with, as user-data.sh wrote it ("" for none)
exit_status() {
  run_log | grep -o "run_all.R exited with status [0-9]*" | tail -1 | grep -o "[0-9]*$"
}

launch() { # prints the new instance id, or nothing
  aws ec2 create-fleet --region "$machine_region" --cli-input-json "file://$fleet" \
    --query 'Instances[0].InstanceIds[0]' --output text 2>/dev/null | grep '^i-'
}

while true; do
  state=$(aws ec2 describe-instances --region "$machine_region" --instance-ids "$id" \
    --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null)
  case "$state" in
    running|pending|"") ;;  # "" = AWS not reachable: try again next time
    shutting-down|stopping) ;;  # the log is still being copied: look next time
    *)
      status=$(exit_status)
      if [ "$status" = "0" ]; then
        echo "$(date +%H:%M) $run finished"
        run_log | tail -20
        exit 0
      fi
      if [ "$status" = "1" ] && run_log | grep -q "Execution halted"; then
        echo "$(date +%H:%M) $run: run_all.R stopped on an error; not launching again"
        run_log | tail -40
        exit 1
      fi
      if [ $relaunches -ge 12 ]; then
        echo "$(date +%H:%M) $run stopped again after 12 new machines: leaving it"
        exit 1
      fi
      echo "$(date +%H:%M) $run: $id is $state without finishing (status: ${status:-none}); launching a new one"
      new=""
      for _ in $(seq 1 20); do  # the spot limit can take a few minutes to free up
        new=$(launch); [ -n "$new" ] && break; sleep 60
      done
      if [ -n "$new" ]; then
        relaunches=$((relaunches + 1)); id=$new
        echo "$(date +%H:%M) $run: new machine $id ($relaunches of at most 12)"
      else
        echo "$(date +%H:%M) $run: no machine after 20 tries; trying again later"
      fi
      ;;
  esac
  sleep "$wait_seconds"
done
