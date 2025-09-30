"""
    default_subdivision(domain::DOM) where {DOM<:AbstractDomain}

Return the default algorithm to subdivide `domain`.
- dimension 1:
    - [`Segment`](@ref): [`subdivide_segment`](@ref)
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
default_subdivision(::Segment) = subdivide_segment
default_subdivision(::Simplex) = subdivide_simplex
default_subdivision(::Triangle) = subdivide_triangle
default_subdivision(::Tetrahedron) = subdivide_tetrahedron
default_subdivision(::Orthotope) = subdivide_orthotope
default_subdivision(::Rectangle) = subdivide_rectangle
default_subdivision(::Cuboid) = subdivide_cuboid

"""
    default_embedded_cubature(domain::DOM) where {DOM<:AbstractDomain}

Return a default embedded cubature for the domains:
- dimension 1:
    - [`Segment`](@ref): [`SEGMENT_GK15`](@ref)
- dimension 2:
    - [`Triangle`](@ref): [`RadonLaurie`](@ref)`()`
    - [`Rectangle`](@ref): [`SQUARE_CH25`](@ref)
- dimension 3:
    - [`Tetrahedron`](@ref): [`GrundmannMoeller`](@ref)`{3}(7, 5)`
    - [`Cuboid`](@ref): [`CUBE_BE65`](@ref)
- dimension `D`:
    - [`Simplex`](@ref): [`GrundmannMoeller`](@ref)`{D}(7, 5)`
    - [`Orthotope`](@ref): [`GenzMalik`](@ref)`{D}()`
"""
@generated function default_embedded_cubature(::Segment{T}) where {T}
    ec = embedded_cubature(SEGMENT_GK15, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    ec = embedded_cubature(GrundmannMoeller{D}(7, 5), T)
    return :($ec)
end
@generated function default_embedded_cubature(::Triangle{T}) where {T}
    ec = embedded_cubature(RadonLaurie(), T)
    return :($ec)
end
@generated function default_embedded_cubature(::Tetrahedron{T}) where {T}
    ec = embedded_cubature(GrundmannMoeller{3}(7, 5), T)
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
