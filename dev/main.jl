using BenchmarkTools
using DataStructures
using LinearAlgebra
using StaticArrays
import HAdaptiveIntegration as hai

abstract type AbstractDomain{D,T} end

struct Orthotope{D,T} <: AbstractDomain{D,T} end

struct Simplex{D,N,T} <: AbstractDomain{D,T} end

const Triangle{T} = Simplex{2,3,T}

abstract type AbstractRule{DOM<:AbstractDomain} end

struct TabulatedRule{DOM<:AbstractDomain} <: AbstractRule{DOM} end

struct GM{D,N,T} <: AbstractRule{Simplex{D,N,T}}
    deg::Int

    function GM{D,T}(deg::Int) where {D,T}
        N = D + 1
        return new{N,D,T}(deg)
    end
end

function fct(domain::DOM, rule::TabulatedRule{DOM}) where {DOM<:AbstractDomain}
    return 0
end

dom = Triangle{Float64}()
rule = TabulatedRule{Triangle}()
fct(dom, rule)
