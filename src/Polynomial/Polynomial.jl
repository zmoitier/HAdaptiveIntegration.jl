module Polynomial

export integral_monomials,
    integral_chebyshev, compute_error_monomials, compute_error_chebyshev

using ..HAdaptiveIntegration.Domain
using ..HAdaptiveIntegration.Rule
using Base.Iterators: countfrom
using Printf: @printf

abstract type AbstractBasis end

struct MonomialBasis <: AbstractBasis end

struct ChebyshevBasis <: AbstractBasis end

include("integrals.jl")
include("errors.jl")

end
