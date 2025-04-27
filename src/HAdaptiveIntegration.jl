module HAdaptiveIntegration

using DataStructures: BinaryHeap
using LinearAlgebra: det, norm
using StaticArrays: SVector, setindex

# This function is intended for internal use to determine the floating-point type.
function promote_to_float(containers...)
    # Determine the type of each container
    container_types = map(eltype, containers)

    # Promote the types of all containers to a common type
    promoted_type = reduce(promote_type, container_types)

    # Convert the promoted type to a floating-point type
    return float(promoted_type)
end

# Supported integration domains
include("Domain/Domain.jl")
using .Domain
export Domain, Simplex, Triangle, Tetrahedron, Orthotope, Rectangle, Cuboid

# Tabulated cubature rule for supported domains
include("Rule/Rule.jl")
using .Rule
export Rule

# Default subdivision and embedded cubature for supported domain
include("default.jl")

# Compute integrals
include("integrate.jl")
export integrate

end
