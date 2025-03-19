using StaticArrays

abstract type Domain{D,T} end

struct Simplex{D,T,N} <: Domain{D,T} end

const Triangle{T} = Simplex{2,T,3}

abstract type Rule{DOM<:Domain} end

struct TabulatedRule{DOM} <: Rule{DOM} end

function fct(domain::DOM, rule::TabulatedRule{DOM}) where {DOM<:Domain} end
