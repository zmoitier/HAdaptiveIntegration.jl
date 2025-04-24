"""
   struct GrundmannMoeller{D} <: AbstractRule{Simplex{D}}

Embedded cubature rule for a `D`-simplex.

## Type Parameters:
- `D`: The dimension of the simplex.

## Fields:
- `order_high::Int`: the high order of the rule.
- `order_low::Int`: the low order of the rule.

## Invariants (check at construction):
- `order_high` and `order_low` must be odd.
- must have `order_high > order_low ≥ 1`.
"""
struct GrundmannMoeller{D} <: AbstractRule{Simplex{D}}
    order_high::Int
    order_low::Int

    function GrundmannMoeller{D}(order_high::Int, order_low::Int) where {D}
        @assert isodd(order_high) && isodd(order_low)
        @assert order_high > order_low ≥ 1
        return new{D}(order_high, order_low)
    end
end

function orders(gm::GrundmannMoeller)
    return gm.order_high, gm.order_low
end
