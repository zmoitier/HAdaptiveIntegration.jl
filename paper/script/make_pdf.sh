#!/usr/bin/bash
set -e

script_dir="$(cd "$(dirname "$0")" && pwd)"
paper_dir="$(cd "$script_dir/.." && pwd)"

for arg in "$@"; do
    case $arg in
        --makefigs)
        echo "=== Generating figures ==="
        bash "$script_dir/make_figures.sh"
        ;;
    esac
done

docker run --rm \
    --volume "$paper_dir":/data \
    --env JOURNAL=joss \
    openjournals/inara

# --user $(id -u):$(id -g) \
