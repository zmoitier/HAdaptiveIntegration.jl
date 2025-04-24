module Domain

"""
    abstract type AbstractDomain{D,T}

Abstract type for integration domains' in `D` dimensions with element type `T`.

## Type Parameters:
- `D`: dimension of the domain.
- `T`: element type of the domain.

## Mandatory methods:
- [`map_from_reference`](@ref)
- [`abs_det_jac`](@ref)

## Useful (but non-mandatory) methods:
- [`reference_domain`](@ref)
- [`map_to_reference`](@ref)
"""
abstract type AbstractDomain{D,T} end

"""
    dimension(::Type{DOM}) where {DOM<:AbstractDomain{D}}

Return the dimension `D` of the given domain `DOM`.
"""
dimension(::Type{AbstractDomain{D,T}}) where {D,T} = D
dimension(::Type{AbstractDomain{D}}) where {D} = D
dimension(::Type{DOM}) where {DOM<:AbstractDomain} = dimension(supertype(DOM))

end
