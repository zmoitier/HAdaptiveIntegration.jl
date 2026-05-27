#!/usr/bin/bash
set -e

for arg in "$@"; do
    case $arg in
        --makefigs)
        echo "=== Generating figures ==="
        bash "$(dirname "$0")/make_figures.sh"
        ;;
    esac
done

docker run --rm \
    --volume $PWD/paper:/data \
    --env JOURNAL=joss \
    openjournals/inara

# --user $(id -u):$(id -g) \
