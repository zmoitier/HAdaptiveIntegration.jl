"""
    struct Simplex{D,N,T} <: AbstractDomain{D,T}

A simplex in `D` dimensions with `N=D+1` vertices of element type `T`.

## Fields:
- `vertices::SVector{N,SVector{D,T}}`: vertices of the simplex.

## Invariants (**not** check at construction):
- `N = D+1`
"""
struct Simplex{D,N,T} <: AbstractDomain{D,T}
    vertices::SVector{N,SVector{D,T}}
end

function Simplex{D,N,T}(vertices::Vararg{SVector{D,T},N}) where {D,N,T}
    return Simplex(SVector{N}(vertices...))
end

"""
    simplex(vertices...)
    simplex(T::DataType, vertices...)

Return a `D`-simplex with element type `T` from a collection of `N=D+1` vertices. Note that
all vertices must have the same length `D`.
"""
function simplex(T::DataType, vertices...)
    @assert allequal(length, vertices) "all `vertices` must have the same length."
    D = length(first(vertices))

    N = D + 1
    @assert length(vertices) == N "Expected $N vertices, but got $(length(vertices))."

    return Simplex(SVector{N}(SVector{D,T}.(vertices)))
end
simplex(vertices...) = simplex(promote_to_float(vertices...), vertices...)

"""
    reference_simplex(D::Int)
    reference_simplex(T::DataType, D::Int)

Return the reference `D`-dimensional simplex with element type `T`, which is the convex hull
of the `N=D+1` points `(0,...,0)`, `(1,0,...,0)`, `(0,1,0,...,0)`, ..., `(0,...,0,1)`.
"""
function reference_simplex(T::DataType, D::Int)
    vertices = [
        zeros(SVector{D,T}), collect(setindex(zeros(SVector{D,T}), 1, i) for i in 1:D)...
    ]
    return Simplex(SVector{D + 1}(vertices))
end

"""
    map_from_reference(domain::DOM) where {DOM<:AbstractDomain}

Return an anonymous function that maps the reference domain to the physical domain `domain`.
"""
function map_from_reference(s::Simplex{D,N,T}) where {D,N,T}
    return u -> begin
        v = (1 - sum(u)) * s.vertices[1]
        for i in 2:N
            v += u[i - 1] * s.vertices[i]
        end
        return v
    end
end

"""
    abs_det_jac(domain::DOM) where {DOM<:AbstractDomain}

Return the absolute value of the Jacobian's determinant of the map from the reference domain
to the physical domain `domain`.
"""
function abs_det_jac(s::Simplex{D,N,T}) where {D,N,T}
    vertices = s.vertices
    jacobian_matrix = hcat(ntuple(i -> vertices[i + 1] - vertices[1], D)...)
    return abs(det(jacobian_matrix))
end

"""
    map_to_reference(domain::DOM) where {DOM<:AbstractDomain}

Return an anonymous function that maps the physical domain `domain` to the reference domain.

## Constraints:
- For `Orthotope`, must have `high_corner .> low_corner`.
- For `Simplex{D}`, the vertices must form a valid `D`-dimensional simplex with non-zero
  volume).
"""
function map_to_reference(s::Simplex{D,N,T}) where {D,N,T}
    v = s.vertices
    jacobian_matrix = hcat(ntuple(i -> v[i + 1] - v[1], D)...)
    @assert !isapprox(det(jacobian_matrix), 0; atol=(√eps(float(T)))) "degenerate $D-dimensional Simplex: the Jacobian matrix is not invertible."

    M = inv(jacobian_matrix)
    return u -> M * (u - v[1])
end

"""
    combinations(n::Int, k::Int)

Helper function to generate all combinations of `k` elements from `1:n`, similar to calling
`combinations(1:n, k)` from `Combinatorics.jl`.
"""
function combinations(n::Int, k::Int)
    function combinations_helper(start::Int, end_::Int, k_::Int)
        if k_ == 0
            return Vector{Int}[[]]
        elseif k_ > (end_ - start + 1)
            return Vector{Int}[]
        else
            res = Vector{Int}[]
            for i in start:(end_ - k_ + 1)
                for c in combinations_helper(i + 1, end_, k_ - 1)
                    push!(res, vcat([i], c))
                end
            end
            return res
        end
    end
    return combinations_helper(1, n, k)
end

"""
    subdivide_reference_simplex(::Val{D}, ::Type{T}=Float64) where {D,T}

Like `subdivide_simplex`, but operates on the reference simplex. Since the output depends
only on the dimension `D`, and the type `T` used to represent coordinates, this function is
generated for each combination of `D` and `T`.
"""
@generated function subdivide_reference_simplex(::Val{D}, (::Type{T})=Float64) where {D,T}
    # vertices of the reference simplex
    vertices = [
        zeros(SVector{D,T}), collect(setindex(zeros(SVector{D,T}), 1, i) for i in 1:D)...
    ]

    # valid permutations of the vertices
    function generate_valid_permutations(n::Int, k::Int)
        valid_perms = []
        # Generate all combinations of k positions from 1:n
        for positions in combinations(n, k)
            # Create a permutation where the first k elements (1:k) are placed at the
            # selected positions in order, and the rest (k+1:n) follow.
            perm = zeros(Int, n)
            selected_pos = sort(positions)
            remaining_pos = setdiff(1:n, selected_pos)
            perm[selected_pos] = 1:k
            perm[remaining_pos] = (k + 1):n
            push!(valid_perms, perm)
        end
        return valid_perms
    end

    sub_simplices = []
    for k in 0:D
        # Compute initial vertex v₀ = (x⁰ + xᵏ) / 2
        x₀ = vertices[1]
        xₖ = vertices[k + 1]
        v₀ = (x₀ .+ xₖ) ./ 2
        # Generate valid permutations for this k
        perms = generate_valid_permutations(D, k)
        for π in perms
            new_vertices = [v₀]
            current_v = v₀
            for ℓ in 1:D
                edge_num = π[ℓ]  # Edge between x^{edge_num-1} and x^{edge_num}
                edge_start = vertices[edge_num]
                edge_end = vertices[edge_num + 1]
                edge_vector = edge_end .- edge_start
                current_v = current_v .+ 0.5 .* edge_vector
                push!(new_vertices, current_v)
            end
            push!(sub_simplices, new_vertices)
        end
    end

    # convert to an efficient format with known sizes
    static_sub_simplices = ntuple(2^D) do i
        pts = SVector{D + 1,SVector{D,T}}(sub_simplices[i])
        Simplex(pts)
    end

    return :($static_sub_simplices)
end

"""
    subdivide_simplex(s::Simplex)

Subdivide a `D`-simplex into `2ᴰ` simplices by using the Freudenthal triangulation.

Implements the `RedRefinementND` algorithm in [Simplicial grid refinement: on Freudenthal's
algorithm and the optimal number of congruence
classes](https://link.springer.com/article/10.1007/s002110050475).
"""
function subdivide_simplex(s::Simplex{D,N,T}) where {D,N,T}
    refs = subdivide_reference_simplex(Val(D), T)
    f = map_from_reference(s)
    map(refs) do ref
        Simplex(f.(ref.vertices))
    end
end

"""
    Triangle{T} = Simplex{2,3,T}

Alias for a 2-dimensional triangle with 3 vertices of value type `T`.
"""
const Triangle{T} = Simplex{2,3,T}

"""
    triangle(a, b, c)
    triangle(T::DataType, a, b, c)

Return a triangle in 2 dimensions, with element type `T`, given by three 2d-points `a`, `b`,
and `c`.
"""
function triangle(T::DataType, a, b, c)
    @assert length(a) == length(b) == length(c) == 2 "`a`, `b`, and `c` must be 2d-vector."
    return simplex(T, a, b, c)
end
triangle(a, b, c) = triangle(promote_to_float(a, b, c), a, b, c)

"""
    subdivide_triangle2(s::Triangle)

Divide the triangle `t` into two triangles by connecting the first point of `t` to the
midpoints of the two other points.
"""
function subdivide_triangle2(t::Triangle{T}) where {T}
    a, b, c = t.vertices
    bc = (b + c) / 2
    return (Triangle{T}(bc, a, b), Triangle{T}(bc, c, a))
end

"""
    subdivide_triangle4(t::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the edges.
"""
function subdivide_triangle4(t::Triangle{T}) where {T}
    a, b, c = t.vertices
    ab = (a + b) / 2
    ac = (c + a) / 2
    bc = (b + c) / 2
    return (
        Triangle{T}(a, ab, ac),
        Triangle{T}(b, bc, ab),
        Triangle{T}(c, ac, bc),
        Triangle{T}(ab, bc, ac),
    )
end

"""
    Tetrahedron{T} = Simplex{3,4,T}

Alias for a 3-dimensional tetrahedron with 4 vertices of element type `T`.
"""
const Tetrahedron{T} = Simplex{3,4,T}

"""
    tetrahedron(a, b, c, d)
    tetrahedron(T::DataType, a, b, c, d)

Return a tetrahedron in 3 dimensions, with element type `T`, given by four 3d-points `a`,
`b`, `c`, and `d`.
"""
function tetrahedron(T::DataType, a, b, c, d)
    @assert length(a) == length(b) == length(c) == length(d) == 3 "`a`, `b`, `c`, and `d` must be 3d-vector."
    return simplex(T, a, b, c, d)
end
tetrahedron(a, b, c, d) = tetrahedron(promote_to_float(a, b, c, d), a, b, c, d)

"""
    subdivide_tetrahedron8(t::Tetrahedron)

Divide the tetrahedron `t` into eight tetrahedra by connecting the midpoints of the edges.
"""
function subdivide_tetrahedron8(t::Tetrahedron{T}) where {T}
    a, b, c, d = t.vertices
    ab = (a + b) / 2
    ac = (a + c) / 2
    ad = (a + d) / 2
    bc = (b + c) / 2
    bd = (b + d) / 2
    cd = (c + d) / 2
    return (
        # (1/2)-tetrahedron on each vertices
        Tetrahedron{T}(a, ab, ac, ad),
        Tetrahedron{T}(ab, b, bc, bd),
        Tetrahedron{T}(ac, bc, c, cd),
        Tetrahedron{T}(ad, bd, cd, d),
        # octahedron splitting
        Tetrahedron{T}(ab, ac, ad, bd),
        Tetrahedron{T}(ab, ac, bc, bd),
        Tetrahedron{T}(ac, ad, bd, cd),
        Tetrahedron{T}(ac, bc, bd, cd),
    )
end
