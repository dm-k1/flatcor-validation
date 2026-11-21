#!/bin/sh
## =====================================================================
##  Wait for an image build of this branch that is still running
##
##  When a push changes the environment (cluster/Dockerfile, its ignore
##  file, renv.lock), .github/workflows/image.yml builds the image for it,
##  and that takes a while. Everything that runs inside the image (the
##  static check, the smoke test, the BridgeStan check, the run-time
##  estimate) first calls this script, so it uses the new image, not the
##  one before it or none at all. It waits while the branch's latest image
##  build is queued or running; it fails if that build failed.
##
##  ci/wait_for_image.sh [seconds]   seconds to wait before the first look,
##                                   so an image build started by the same
##                                   push shows up (default 60)
##  Needs GH_TOKEN, GITHUB_REPOSITORY and GITHUB_REF_NAME (set in Actions).
## =====================================================================
set -e
sleep "${1:-60}"
while true; do
  state=$(gh run list -R "$GITHUB_REPOSITORY" --workflow image.yml \
    --branch "$GITHUB_REF_NAME" --limit 1 --json status,conclusion \
    -q '.[0] // empty | "\(.status) \(.conclusion)"')
  case "$state" in
    "") echo "No image build for this branch"; exit 0 ;;
    "completed success") echo "Image built"; exit 0 ;;
    completed*) echo "Image build ended: $state"; exit 1 ;;
  esac
  sleep 60
done
