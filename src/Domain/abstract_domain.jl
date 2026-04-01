"""
    abstract type AbstractDomain{D,T}

Abstract type for integration domains in `D` dimensions with element type `T`.

## Type Parameters:
- `D`: dimension of the domain.
- `T`: element type of the domain.

## Mandatory methods:
- [`map_from_reference`](@ref)

## Useful (but non-mandatory) methods:
- [`reference_domain`](@ref)
- [`map_to_reference`](@ref)
"""
abstract type AbstractDomain{D, T} end

"""
    Φ, μ = map_from_reference(domain::DOM) where {DOM<:AbstractDomain}

Return `(Φ, μ)`, where `Φ` maps the reference domain to the physical domain `domain`, and
`μ` is the absolute value of the Jacobian determinant of `Φ`.
"""
function map_from_reference end

"""
    map_to_reference(domain::DOM) where {DOM<:AbstractDomain}

Return an anonymous function that maps the physical `domain` to the reference domain.

## Constraints:
- For `Segment`, must have `xmax > xmin`.
- For `Simplex{D}`, the vertices must form a valid `D`-dimensional simplex with non-zero
  volume.
- For `Orthotope`, must have `high_corner .> low_corner`.
"""
function map_to_reference end

"""
    reference_domain(::Type{<:AbstractDomain})

Return the reference domain for the given domain type. See [`reference_segment`](@ref),
[`reference_simplex`](@ref), and [`reference_orthotope`](@ref) for more details.
"""
function reference_domain end
