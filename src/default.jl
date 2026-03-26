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
- dimension `D`:
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
function default_embedded_cubature(::Segment{T}) where {T}
    S = typeof(one(T))
    return cache_dec(Segment{S}, SEGMENT_GK15)
end
function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    S = typeof(one(T))
    if D == 2
        return cache_dec(Triangle{S}, RadonLaurie())
    end
    return cache_dec(Simplex{D,S,N}, GrundmannMoeller{D}(7, 5))
end
function default_embedded_cubature(::Orthotope{D,T}) where {D,T}
    S = typeof(one(T))
    if D == 2
        return cache_dec(Rectangle{S}, SQUARE_CH25)
    end
    if D == 3
        return cache_dec(Cuboid{S}, CUBE_BE65)
    end
    return cache_dec(Orthotope{D,S}, GenzMalik{D}())
end

# cache for default embedded cubatures
const CACHE_DEC = Dict{DataType,EmbeddedCubature}()
const CACHE_DEC_LOCK = ReentrantLock()

function cache_dec(
    ::Type{DOM}, rule::AR
)::EmbeddedCubature{D,T} where {D,T,DOM<:AbstractDomain{D,T},AR<:AbstractRule}
    lock(CACHE_DEC_LOCK)
    try
        return get!(CACHE_DEC, DOM) do
            embedded_cubature(rule, T)
        end
    finally
        unlock(CACHE_DEC_LOCK)
    end
end
