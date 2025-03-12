"""
    default_subdivision(domain::Domain)

Return the default algorithm to subdivide `domain`.
- dimension 1:
    - [`segment`](@ref): [`subdivide_segment2`](@ref)
- dimension 2:
    - [`rectangle`](@ref): [`subdivide_rectangle4`](@ref)
    - [`triangle`](@ref): [`subdivide_triangle4`](@ref)
- dimension 3:
    - [`cuboid`](@ref): [`subdivide_cuboid8`](@ref)
    - [`tetrahedron`](@ref): [`subdivide_tetrahedron8`](@ref)
- dimension `d`:
    - [`simplex`](@ref): [`subdivide_simplex`](@ref)
"""
default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Rectangle) = subdivide_rectangle4
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Cuboid) = subdivide_cuboid8
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Simplex) = subdivide_simplex

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
    ec = embedded_cubature(T, SEGMENT_GK15)
    return :($ec)
end
@generated function default_embedded_cubature(::Rectangle{T}) where {T}
    ec = embedded_cubature(T, SQUARE_CHG25)
    return :($ec)
end
@generated function default_embedded_cubature(::Triangle{T}) where {T}
    ec = embedded_cubature(T, TRIANGLE_RL19)
    return :($ec)
end
@generated function default_embedded_cubature(::Cuboid{T}) where {T}
    ec = embedded_cubature(T, CUBE_BE65)
    return :($ec)
end
@generated function default_embedded_cubature(::Tetrahedron{T}) where {T}
    ec = embedded_cubature(T, TETRAHEDRON_GM35)
    return :($ec)
end
@generated function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    ec = embedded_cubature(T, GrundmannMoeller(D, 7))
    return :($ec)
end
