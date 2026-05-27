#!/usr/bin/bash
set -e

dir="$(cd "$(dirname "$0")" && pwd)"

for script in cvg_triangle.jl cvg_tetrahedron.jl cvg_rectangle.jl cvg_cuboid.jl; do
    echo ""
    echo "============================================================"
    echo "Generating figures: $script"
    echo "============================================================"
    julia --project="$dir" "$dir/$script"
done

echo ""
echo "All figures generated."
