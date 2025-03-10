module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

using GrundmannMoeller: grundmann_moeller

# Supported integration domains
include("domain.jl")
export orthotope, segment, rectangle, cuboid, simplex, triangle, tetrahedron

# Subdivision strategies for various domains
include("domain_subdivision.jl")

# Embedded cubature
include("cubature_embedded.jl")
include("cubature_check.jl")

# Tabulated cubature rule for supported domains
include("rule_orthotope.jl")
include("rule_simplex.jl")

"""
    default_embedded_cubature(domain::Domain)

Return a default embedded cubature for the domains:
- dimension 1:
    - [`segment`](@ref): [`SEGMENT_GK15`](@ref)
- dimension 2:
    - [`rectangle`](@ref): [`SQUARE_CHG25`](@ref)
    - [`triangle`](@ref): [`TRIANGLE_RL19`](@ref)
- dimension 3:
    - [`cuboid`](@ref): [`CUBE_BE65`](@ref)
    - [`tetrahedron`](@ref): [`TETRAHEDRON_GM35`](@ref)
- dimension `d`:
    - [`simplex`](@ref): `[GrundmannMoeller](@ref)(d, 7)`
"""
@generated function default_embedded_cubature(::Segment{T}) where {T}
    ec = embedded_cubature(SEGMENT_GK15, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Rectangle{T}) where {T}
    ec = embedded_cubature(SQUARE_CHG25, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Triangle{T}) where {T}
    ec = embedded_cubature(TRIANGLE_RL19, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Cuboid{T}) where {T}
    ec = embedded_cubature(CUBE_BE65, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Tetrahedron{T}) where {T}
    ec = embedded_cubature(TETRAHEDRON_GM35, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    ec = embedded_cubature(GrundmannMoeller(D, 7), T)
    return :($ec)
end

# Compute integrals
include("integrate.jl")
export integrate

end
