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

const _DEFAULT_EC_CACHE = Dict{DataType,Any}()
const _DEFAULT_EC_CACHE_LOCK = ReentrantLock()

function _cached_default_ec(
    ::Type{DOM}, rule::AR, (::Type{T})
)::EmbeddedCubature{D,T} where {D,DOM<:AbstractDomain{D},AR<:AbstractRule,T<:Real}
    lock(_DEFAULT_EC_CACHE_LOCK)
    try
        return get!(_DEFAULT_EC_CACHE, DOM) do
            embedded_cubature(rule, T)
        end
    finally
        unlock(_DEFAULT_EC_CACHE_LOCK)
    end
end

_scalar_type(::Type{T}) where {T} = typeof(one(T))

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
    S = _scalar_type(T)
    return _cached_default_ec(Segment{T}, SEGMENT_GK15, S)
end
function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    S = _scalar_type(T)
    return _cached_default_ec(Simplex{D,T,N}, GrundmannMoeller{D}(7, 5), S)
end
function default_embedded_cubature(::Triangle{T}) where {T}
    S = _scalar_type(T)
    return _cached_default_ec(Triangle{T}, RadonLaurie(), S)
end
function default_embedded_cubature(::Tetrahedron{T}) where {T}
    S = _scalar_type(T)
    return _cached_default_ec(Tetrahedron{T}, GrundmannMoeller{3}(7, 5), S)
end
function default_embedded_cubature(::Orthotope{D,T}) where {D,T}
    S = _scalar_type(T)
    return _cached_default_ec(Orthotope{D,T}, GenzMalik{D}(), S)
end
function default_embedded_cubature(::Rectangle{T}) where {T}
    S = _scalar_type(T)
    return _cached_default_ec(Rectangle{T}, SQUARE_CH25, S)
end
function default_embedded_cubature(::Cuboid{T}) where {T}
    S = _scalar_type(T)
    return _cached_default_ec(Cuboid{T}, CUBE_BE65, S)
end
