# TODO: move to utils
"""
    check_subdivision(
        domain::Domain{D,T},
        subdiv_algo;
        atol::T=zero(T),
        rtol::T=(atol > zero(T)) ? zero(T) : 10 * eps(T),
    ) where {D,T<:Real}

Return `nothing` if the sum of the volume of the subdomain by the `subdiv_algo` is equal to the volume of the domain else throw an error.
"""
function check_subdivision(
    domain::Domain{D,T},
    subdiv_algo;
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T<:Real}
    subdomains = subdiv_algo(domain)
    if isapprox(
        sum(abs_det_jacobian.(subdomains)), abs_det_jacobian(domain); atol=atol, rtol=rtol
    )
        @info "`$(Symbol(subdiv_algo))` pass volume test."
        return nothing
    else
        @error "`$(Symbol(subdiv_algo))` do not partition the domain."
    end
end

"""
    subdivide_segment2(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment2(s::Segment{T}) where {T<:Real}
    a, b = s.low_corner, s.high_corner
    m = (a + b) / 2
    return (Segment{T}(a, m), Segment{T}(m, b))
end

"""
    subdivide_segment3(s::Segment)

Divide the segment `s` into three segments of equal length.
"""
function subdivide_segment3(s::Segment{T}) where {T<:Real}
    a, b = s.low_corner, s.high_corner
    m1, m2 = (2 * a + b) / 3, (a + 2 * b) / 3
    return (Segment{T}(a, m1), Segment{T}(m1, m2), Segment{T}(m2, b))
end

"""
    subdivide_triangle2(s::Triangle)

Divide the triangle `t` into two triangles by connecting the first point of `t` to the midpoints of
the two other points.
"""
function subdivide_triangle2(t::Triangle{T}) where {T<:Real}
    a, b, c = t.points
    bc = (b + c) / 2
    return (Triangle{T}(bc, a, b), Triangle{T}(bc, c, a))
end

"""
    subdivide_triangle4(t::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the edges.
"""
function subdivide_triangle4(t::Triangle{T}) where {T<:Real}
    a, b, c = t.points
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
    subdivide_rectangle4(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the midpoints
of the edges.
"""
function subdivide_rectangle4(r::Rectangle{T}) where {T<:Real}
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
    subdivide_tetrahedron8(t::Tetrahedron)

Divide the tetrahedron `t` into eight tetrahedra by connecting the midpoints of the edges.
"""
function subdivide_tetrahedron8(t::Tetrahedron{T}) where {T<:Real}
    a, b, c, d = t.points
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
        Tetrahedron{T}(ab, bc, bd, ad),
        Tetrahedron{T}(ac, bc, cd, ad),
        Tetrahedron{T}(ad, cd, bd, bc),
        Tetrahedron{T}(bc, ac, ab, ad),
    )
end

"""
    _subdivide_reference_simplex_freudenthal(::Val{D}, ::Type{T}=Float64)

Like `subdivide_simplex_freudenthal`, but operates on the reference simplex. Since the
output depends only on the dimension `D`, and the type `T` used to represent coordinates,
this function is generated for each combination of `D` and `T`.
"""
@generated function _subdivide_reference_simplex_freudenthal(
    ::Val{D}, ::Type{T}=Float64
) where {D,T}
    # vertices of the reference simplex
    vertices = [zeros(T, D)]
    for i in 1:D
        v = zeros(T, D)
        v[i] = 1
        push!(vertices, v)
    end

    # valid permutations of the vertices
    function generate_valid_permutations(n, k)
        valid_perms = []
        # Generate all combinations of k positions from 1:n
        for positions in combinations(1:n, k)
            # Create a permutation where the first k elements (1:k) are placed
            # at the selected positions in order, and the rest (k+1:n) follow.
            perm = zeros(Int, n)
            selected_pos = sort(positions)
            remaining_pos = setdiff(1:n, selected_pos)
            perm[selected_pos] = 1:k
            perm[remaining_pos] = (k + 1):n
            push!(valid_perms, perm)
        end
        return valid_perms
    end

    n = length(vertices) - 1  # T is an n-simplex with (n+1) vertices
    subsimplices = []
    for k in 0:n
        # Compute initial vertex v₀ = (x⁰ + xᵏ) / 2
        x₀ = vertices[1]
        xₖ = vertices[k + 1]
        v₀ = (x₀ .+ xₖ) ./ 2
        # Generate valid permutations for this k
        perms = generate_valid_permutations(n, k)
        for π in perms
            new_vertices = [v₀]
            current_v = v₀
            for ℓ in 1:n
                edge_num = π[ℓ]  # Edge between x^{edge_num-1} and x^{edge_num}
                edge_start = vertices[edge_num]
                edge_end = vertices[edge_num + 1]
                edge_vector = edge_end .- edge_start
                current_v = current_v .+ 0.5 .* edge_vector
                push!(new_vertices, current_v)
            end
            push!(subsimplices, new_vertices)
        end
    end
    # convert to an efficient format with known sizes

    static_subsimplices = ntuple(2^D) do i
        pts = SVector{D + 1,SVector{D,T}}(subsimplices[i])
        Simplex(pts)
    end

    return :($static_subsimplices)
end

"""
    subdivide_simplex_freudenthal(s::Simplex)

Subdivive a `D`-dimensional simplex into `2ᴰ` simplices by using the Freudenthal
triangulation.

Implements the `RedRefinementND`` algorithm in [Simplicial grid refinement: on Freudenthal's
algorithm and the optimal number of congruence
classes](https://link.springer.com/article/10.1007/s002110050475).
"""
function subdivide_simplex_freudenthal(s::Simplex{D,T}) where {D,T}
    refs = _subdivide_reference_simplex_freudenthal(Val(D), T)
    f = map_from_reference(s)
    map(refs) do ref
        Simplex(f.(ref.points))
    end
end

"""
    subdivide_cuboid8(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints of the
edges.
"""
function subdivide_cuboid8(c::Cuboid{T}) where {T<:Real}
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

const LIST_SUBDIVISION_ALGO = [
    "segment" => ["subdivide_segment2", "subdivide_segment3"],
    "rectangle" => ["subdivide_rectangle4"],
    "triangle" => ["subdivide_triangle2", "subdivide_triangle4"],
    "cuboid" => ["subdivide_cuboid8"],
    "tetrahedron" => ["subdivide_tetrahedron8"],
]

function default_subdivision(d::Domain)
    @error "no default subdivision for $(typeof(d))."
end
default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Rectangle) = subdivide_rectangle4
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Cuboid) = subdivide_cuboid8
