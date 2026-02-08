using BenchmarkTools
using HAdaptiveIntegration:
    Simplex,
    Tetrahedron,
    Triangle,
    subdivide_simplex,
    subdivide_tetrahedron,
    subdivide_triangle,
    Domain
using StaticArrays

function subdivide_simplex_scheme(::Val{D}) where {D}
    N = D + 1
    K = 1 << D # 2^D
    color_schemes = MVector{K,SMatrix{2,N,Int}}(undef)

    color_scheme = MMatrix{2,N,Int}(undef)
    for k in 0:(K - 1)
        color = 1

        color_scheme[1, 1] = color
        for j in 1:D
            if ((k >> (j - 1)) & 1) == 0
                color += 1
            end
            color_scheme[1, j + 1] = color
        end

        color_scheme[2, 1] = color
        for j in 1:D
            if ((k >> (j - 1)) & 1) == 1
                color += 1
            end
            color_scheme[2, j + 1] = color
        end

        color_schemes[k + 1] = SMatrix(color_scheme)
    end

    return SVector(color_schemes)
end

# function subdivide_simplex_scheme2(d::Int)
#     N = 1 << d # 2^d
#     color_schemes = Vector{SMatrix{2,d + 1,Int}}(undef, N)

#     color_scheme = zeros(MMatrix{2,d + 1,Int})
#     for n in 0:(N - 1)
#         c1 = 1
#         c2 = 0

#         color_scheme[1, 1] = c1
#         color_scheme[2, 1] = c2
#         for j in 1:d
#             if ((n >> (j - 1)) & 1) == 0
#                 c1 += 1
#             else
#                 c2 += 1
#             end
#             color_scheme[1, j + 1] = c1
#             color_scheme[2, j + 1] = c2
#         end
#         color_scheme[2, :] .+= c1

#         color_schemes[n + 1] = color_scheme
#     end

#     return color_schemes
# end

function main()
    # triangle = Triangle((0.10, 0.44), (0.01, 0.29), (0.94, 0.99))
    # display(@benchmark subdivide_simplex($triangle))
    # display(@benchmark subdivide_triangle($triangle))

    # tetrahedron = Tetrahedron(
    #     (0.10, 0.44, 0.27), (0.01, 0.29, 0.67), (0.94, 0.99, 0.81), (0.15, 0.65, 0.57)
    # )
    # display(@benchmark subdivide_simplex($tetrahedron))
    # display(@benchmark subdivide_tetrahedron($tetrahedron))

    @btime Domain.subdivide_reference_simplex(Val(8))

    return nothing
end

main()
