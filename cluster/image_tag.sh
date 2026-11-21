#!/bin/sh
## =====================================================================
##  The tag of the container image for this checkout
##
##  The image holds only the environment (R, the packages of renv.lock,
##  CmdStan, BridgeStan); it is built from exactly three files. Their
##  fingerprint names the image: .github/workflows/image.yml tags every image
##  it publishes with it, and whatever runs the validation (the AWS user
##  data, the CI workflows, a cluster) pulls the image with the tag of its own
##  checkout. So any commit, on any branch, gets the image built for its
##  environment, and a commit whose environment has no image yet fails at the
##  pull instead of running in the wrong one.
##
##  cluster/image_tag.sh  prints  env-<12 hex>  (from any folder)
##  Image: ghcr.io/dm-k1/flatcor-validation:$(cluster/image_tag.sh)
## =====================================================================
set -e
cd "$(dirname "$0")/.."
if command -v sha256sum >/dev/null 2>&1; then
  hash=$(cat cluster/Dockerfile cluster/Dockerfile.dockerignore renv.lock | sha256sum)
else
  hash=$(cat cluster/Dockerfile cluster/Dockerfile.dockerignore renv.lock | shasum -a 256)
fi
echo "env-$(echo "$hash" | cut -c1-12)"
