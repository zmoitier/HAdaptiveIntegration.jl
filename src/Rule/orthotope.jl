"""
   struct GenzMalik{D} <: AbstractRule{Orthotope{D}}

Embedded cubature rule for a `D`-orthotope of high order `7` and low order `5`.

## Type Parameters:
- `D`: The dimension of the orthotope.
"""
struct GenzMalik{D} <: AbstractRule{Orthotope{D}} end

function orders(::GenzMalik)
    return 7, 5
end
