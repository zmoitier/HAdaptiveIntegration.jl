using BenchmarkTools
using DataStructures
using LinearAlgebra
using StaticArrays
import HAdaptiveIntegration as hai

# abstract type AbstractDomain{D} end

# struct Orthotope{D,T} <: AbstractDomain{D} end

# struct Simplex{D,N,T} <: AbstractDomain{D} end

# const Triangle{T} = Simplex{2,3,T}

# abstract type AbstractRule{DOM<:AbstractDomain} end

# struct TabulatedRule{DOM<:AbstractDomain} <: AbstractRule{DOM} end

# struct GM{D,N,T} <: AbstractRule{Simplex{D,N,T}}
#     deg::Int

#     function GM{D,T}(deg::Int) where {D,T}
#         N = D + 1
#         return new{D,N,T}(deg)
#     end
# end

# function fct(domain::DOM, rule::TabulatedRule{DOM}) where {DOM<:AbstractDomain}
#     return 0
# end

# dom = Triangle{Float64}()
# rule = TabulatedRule{Triangle}()
# fct(dom, rule)

fct = x -> 1 / norm(x)
domain = hai.reference_simplex(3)

rtol = √eps(Float64)
I, E = hai.integrate(fct, domain; rtol=rtol)
R = 0.3614258523411 # reference computed using BigFloat
@show abs(I - R) ≤ rtol * abs(R)

display(@benchmark hai.integrate($fct, $domain; rtol=$rtol))
