"""
    check_subdivision(
        subdiv_algo,
        domain::AbstractDomain{D,T};
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
    ) where {D,T}

Return 0 if the sum of the volume of the subdomain by the `subdiv_algo` is equal to the
volume of the domain else return 1.
"""
function check_subdivision(
    subdiv_algo,
    domain::AbstractDomain{D,T};
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T}
    sub_domains = subdiv_algo(domain)
    if !isapprox(sum(abs_det_jac.(sub_domains)), abs_det_jac(domain); atol=atol, rtol=rtol)
        @error "`$(Symbol(subdiv_algo))` do not partition the domain within the tolerance."
        return 1
    end

    @info "`$(Symbol(subdiv_algo))` pass volume test."
    return 0
end

"""
    subdivide_segment2(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment2(s::Segment{T}) where {T}
    a, b = s.low_corner, s.high_corner
    m = (a + b) / 2
    return (Segment{T}(a, m), Segment{T}(m, b))
end

"""
    subdivide_rectangle4(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the
midpoints of the edges.
"""
function subdivide_rectangle4(r::Rectangle{T}) where {T}
    a, b = r.low_corner, r.high_corner
    m = (a + b) / 2
    return (
        Rectangle{T}(a, m),
        Rectangle{T}(SVector(m[1], a[2]), SVector(b[1], m[2])),
        Rectangle{T}(SVector(a[1], m[2]), SVector(m[1], b[2])),
        Rectangle{T}(m, b),
    )
end

"""
    subdivide_cuboid8(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints
of the edges.
"""
function subdivide_cuboid8(c::Cuboid{T}) where {T}
    a, b = c.low_corner, c.high_corner
    m = (a + b) / 2
    return (
        Cuboid{T}(a, m),
        Cuboid{T}(SVector(m[1], a[2], a[3]), SVector(b[1], m[2], m[3])),
        Cuboid{T}(SVector(a[1], m[2], a[3]), SVector(m[1], b[2], m[3])),
        Cuboid{T}(SVector(a[1], a[2], m[3]), SVector(m[1], m[2], b[3])),
        Cuboid{T}(SVector(a[1], m[2], m[3]), SVector(m[1], b[2], b[3])),
        Cuboid{T}(SVector(m[1], a[2], m[3]), SVector(b[1], m[2], b[3])),
        Cuboid{T}(SVector(m[1], m[2], a[3]), SVector(b[1], b[2], m[3])),
        Cuboid{T}(m, b),
    )
end

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
@generated function subdivide_reference_simplex(::Val{D}, ::Type{T}=Float64) where {D,T}
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
