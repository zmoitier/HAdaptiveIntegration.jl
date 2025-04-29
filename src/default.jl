"""
    default_subdivision(domain::DOM) where {DOM<:AbstractDomain}

Return the default algorithm to subdivide `domain`.
- dimension 2:
    - [`Triangle`](@ref): [`subdivide_triangle`](@ref)
    - [`Rectangle`](@ref): [`subdivide_rectangle`](@ref)
- dimension 3:
    - [`Tetrahedron`](@ref): [`subdivide_tetrahedron`](@ref)
    - [`Cuboid`](@ref): [`subdivide_cuboid`](@ref)
- dimension `d`:
    - [`Simplex`](@ref): [`subdivide_simplex`](@ref)
    - [`Orthotope`](@ref): [`subdivide_orthotope`](@ref)
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
    - [`Triangle`](@ref): [`TRIANGLE_RL19`](@ref)
    - [`Rectangle`](@ref): [`SQUARE_CH25`](@ref)
- dimension 3:
    - [`Tetrahedron`](@ref): [`TETRAHEDRON_GM35`](@ref)
    - [`Cuboid`](@ref): [`CUBE_BE65`](@ref)
- dimension `d`:
    - [`Simplex`](@ref): [`GrundmannMoeller{d}(7, 5)`](@ref)`
    - [`Orthotope`](@ref): [`GenzMalik{d}()`](@ref)
"""
@generated function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    ec = embedded_cubature(GrundmannMoeller{D}(7, 5), T)
    return :($ec)
end
@generated function default_embedded_cubature(::Triangle{T}) where {T}
    ec = embedded_cubature(TRIANGLE_RL19, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Tetrahedron{T}) where {T}
    ec = embedded_cubature(TETRAHEDRON_GM35, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Orthotope{D,T}) where {D,T}
    ec = embedded_cubature(GenzMalik{D}(), T)
    return :($ec)
end
@generated function default_embedded_cubature(::Rectangle{T}) where {T}
    ec = embedded_cubature(SQUARE_CH25, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Cuboid{T}) where {T}
    ec = embedded_cubature(CUBE_BE65, T)
    return :($ec)
end
