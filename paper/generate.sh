#!/usr/bin/env bash
set -e

makefigs=false

for arg in "$@"; do
    case $arg in
        makefigs=true)  makefigs=true ;;
        makefigs=false) makefigs=false ;;
    esac
done

if [ "$makefigs" = true ]; then
    echo "=== Generating figures ==="
    bash "$(dirname "$0")/generate_figures.sh"
fi

docker run --rm \
    --volume $PWD/paper:/data \
    --env JOURNAL=joss \
    openjournals/inara

# --user $(id -u):$(id -g) \
