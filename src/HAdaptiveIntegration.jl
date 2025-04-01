module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

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
include("domain.jl")
export orthotope, segment, rectangle, cuboid, simplex, triangle, tetrahedron

# Subdivision strategies for various domains
include("domain_subdivision.jl")

# Tabulated cubature rule for supported domains
include("rule.jl")
include("rule_orthotope.jl")
include("rule_simplex.jl")

# Embedded cubature
include("cubature_embedded.jl")
include("cubature_check.jl")

# Default subdivision and embedded cubature for supported domain
include("default.jl")

# Compute integrals
include("integrate.jl")
export integrate

end
