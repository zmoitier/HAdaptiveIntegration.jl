"""
    abstract type Domain{D,T<:Real}

Abstract type for integration domains in `D` dimensions with type `T<:Real`.
"""
abstract type Domain{D,T<:Real} end

"""
    reference_domain(domain_type::DataType)

Return the reference domain for `domain_type`.
"""
function reference_domain(domain_type::DataType)
    @assert (domain_type <: Domain) "only subtype of `Domain` are allowed."
    D = domain_type.parameters[1]
    T = domain_type.parameters[2]

    if domain_type <: Simplex
        return reference_simplex(D, T)
    end

    if domain_type <: OrthotopeTBW
        return reference_orthotope(D, T)
    end

    @assert false "no reference domain implemented for $domain_type."
end

"""
    map_from_reference(d::Domain)::Function

Return an anonymous function that maps the reference domain to the physical domain `d`.
"""
function map_from_reference(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_from_reference` is not implemented for type $(typeof(d))."
end

"""
    map_to_reference(d::Domain)::Function

Return an anonymous function that maps the physical domain `d` to the reference domain.
"""
function map_to_reference(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_to_reference` is not implemented for type $(typeof(d))."
end

"""
    abs_det_jacobian(d::Domain)

The absolute value of the Jacobian's determinant of the map from the reference domain to the physical domain `d`.
"""
function abs_det_jacobian(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_from_reference` is not implemented for $(typeof(d))."
end
