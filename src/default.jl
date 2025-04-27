"""
    default_subdivision(domain::DOM) where {DOM<:AbstractDomain}

Return the default algorithm to subdivide `domain`.
- dimension 2:
    - [`triangle`](@ref): [`subdivide_triangle`](@ref)
    - [`rectangle`](@ref): [`subdivide_rectangle`](@ref)
- dimension 3:
    - [`tetrahedron`](@ref): [`subdivide_tetrahedron`](@ref)
    - [`cuboid`](@ref): [`subdivide_cuboid`](@ref)
- dimension `d`:
    - [`simplex`](@ref): [`subdivide_simplex`](@ref)
    - [`orthotope`](@ref): [`subdivide_orthotope`](@ref)
"""
default_subdivision(::Simplex) = subdivide_simplex
default_subdivision(::Triangle) = subdivide_triangle
default_subdivision(::Tetrahedron) = subdivide_tetrahedron
default_subdivision(::Orthotope) = subdivide_orthotope
default_subdivision(::Rectangle) = subdivide_rectangle
default_subdivision(::Cuboid) = subdivide_cuboid

"""
    default_embedded_cubature(domain::DOM) where {DOM<:AbstractDomain}

Return a default embedded cubature for the domains:
- dimension 2:
    - [`triangle`](@ref): [`TRIANGLE_RL19`](@ref)
    - [`rectangle`](@ref): [`SQUARE_CH25`](@ref)
- dimension 3:
    - [`tetrahedron`](@ref): [`TETRAHEDRON_GM35`](@ref)
    - [`cuboid`](@ref): [`CUBE_BE65`](@ref)
- dimension `d`:
    - [`simplex`](@ref): [`GrundmannMoeller{d}(7, 5)`](@ref)`
    - [`orthotope`](@ref): [`GenzMalik{d}()`](@ref)
"""
@generated function default_embedded_cubature(::Triangle{T}) where {T}
    ec = embedded_cubature(T, TRIANGLE_RL19)
    return :($ec)
end
@generated function default_embedded_cubature(::Rectangle{T}) where {T}
    ec = embedded_cubature(T, SQUARE_CH25)
    return :($ec)
end
@generated function default_embedded_cubature(::Tetrahedron{T}) where {T}
    ec = embedded_cubature(T, TETRAHEDRON_GM35)
    return :($ec)
end
@generated function default_embedded_cubature(::Cuboid{T}) where {T}
    ec = embedded_cubature(T, CUBE_BE65)
    return :($ec)
end
@generated function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    ec = embedded_cubature(T, GrundmannMoeller{D}(7, 5))
    return :($ec)
end
@generated function default_embedded_cubature(::Orthotope{D,T}) where {D,T}
    ec = embedded_cubature(T, GenzMalik{D}())
    return :($ec)
end
