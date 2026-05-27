#!/usr/bin/env bash
set -e

dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for script in cvg_triangle.jl cvg_tetra.jl cvg_rectangle.jl cvg_cube.jl; do
    echo ""
    echo "============================================================"
    echo "Generating figures: $script"
    echo "============================================================"
    julia --project="$dir" "$dir/$script"
done

echo ""
echo "All figures generated."
