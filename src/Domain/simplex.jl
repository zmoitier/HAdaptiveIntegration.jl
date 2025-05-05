"""
    struct Simplex{D,T,N} <: AbstractDomain{D,T}

A simplex in `D` dimensions with `N=D+1` vertices of element type `T`.

## Fields:
- `vertices::SVector{N,SVector{D,T}}`: vertices of the simplex.

## Invariants (**not** check at construction):
- `N = D+1`

## Constructors:
- `Simplex(vertices...)`
- `Simplex{T}(vertices...)`
- `Simplex(vertices::SVector{N,SVector{D,T}})`
"""
struct Simplex{D,T,N} <: AbstractDomain{D,T}
    vertices::SVector{N,SVector{D,T}}
end

function Simplex{T}(vertices...) where {T}
    N = length(vertices)
    D = N - 1

    for vertex in vertices
        @assert length(vertex) == D "$vertex must have length $D."
    end

    return Simplex(SVector{N}(SVector{D,T}.(vertices)...))
end

function Simplex(vertices...)
    N = length(vertices)
    D = N - 1

    for vertex in vertices
        @assert length(vertex) == D "$vertex must have length $D."
    end

    return Simplex(SVector{N}(float.(SVector{D}.(vertices))...))
end

"""
    reference_simplex(D::Int, T=float(Int))

Return the reference `D`-dimensional simplex with element type `T`, which is the convex hull
of the `N=D+1` points `(0,...,0)`, `(1,0,...,0)`, `(0,1,0,...,0)`, ..., `(0,...,0,1)`.
"""
function reference_simplex(D::Int, (::Type{T})=float(Int)) where {T}
    vertices = [
        zeros(SVector{D,T}), collect(setindex(zeros(SVector{D,T}), 1, i) for i in 1:D)...
    ]
    return Simplex(SVector{D + 1}(vertices))
end

"""
    map_from_reference(domain::DOM) where {DOM<:AbstractDomain}

Return an anonymous function that maps the reference domain to the physical domain `domain`.
"""
function map_from_reference(s::Simplex{D,T,N}) where {D,T,N}
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
function abs_det_jac(s::Simplex{D,T,N}) where {D,T,N}
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
function map_to_reference(s::Simplex{D,T,N}) where {D,T,N}
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
    subdivide_reference_simplex(::Val{D}, ::Type{T}=float(Int)) where {D,T}

Like `subdivide_simplex`, but operates on the reference simplex. Since the output depends
only on the dimension `D`, and the type `T` used to represent coordinates, this function is
generated for each combination of `D` and `T`.
"""
@generated function subdivide_reference_simplex(
    ::Val{D}, (::Type{T})=float(Int)
) where {D,T}
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

    N = D + 1
    sub_simplices = Vector{SVector{N,SVector{D,T}}}()
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
            push!(sub_simplices, SVector{N}(new_vertices))
        end
    end

    # convert to an efficient format with known sizes
    static_sub_simplices = ntuple(2^D) do i
        Simplex(sub_simplices[i])
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
function subdivide_simplex(s::Simplex{D,T,N}) where {D,T,N}
    refs = subdivide_reference_simplex(Val(D), T)
    f = map_from_reference(s)
    map(refs) do ref
        Simplex(f.(ref.vertices))
    end
end
