module Polynomial

export integral_monomials_exact, integral_monomials_rule

using ..HAdaptiveIntegration.Domain
using ..HAdaptiveIntegration.Rule
using Base.Iterators: countfrom
using Printf: @printf

include("integral_exact.jl")
include("integral_num.jl")

end
