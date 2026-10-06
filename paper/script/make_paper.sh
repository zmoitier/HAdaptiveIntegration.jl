#!/usr/bin/bash
set -e

paper_dir="$(cd "$(dirname "$0")/.." && pwd)"

docker run --rm \
    --volume "$paper_dir":/data \
    --env JOURNAL=joss \
    openjournals/inara
