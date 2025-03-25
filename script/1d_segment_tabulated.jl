import HAdaptiveIntegration as hai

using Printf: Format, format
using StaticArrays

struct AffineMap{D,T,N}
    A::SMatrix{D,D,T,N}
    b::SVector{D,T}
end

function (Φ::AffineMap{D,T,N})(x::SVector{D,T}) where {D,T,N}
    return Φ.A * x + Φ.b
end

function assemble(points...)
    matrices = Vector{Vector{SMatrix{1,1,BigFloat,1}}}()
    nodes = Vector{SVector{1,BigFloat}}()
    weights_high = Vector{BigFloat}()
    weights_low = Vector{BigFloat}()

    for (orbit_mat, point, wh, wl) in points
        push!(matrices, orbit_mat)
        push!(nodes, SVector{1}(parse(BigFloat, point)))
        push!(weights_high, parse(BigFloat, wh))
        if wl != ""
            push!(weights_low, parse(BigFloat, wl))
        end
    end

    return (
        matrices=matrices, nodes=nodes, weights_high=weights_high, weights_low=weights_low
    )
end

# [ nodes[1]..., ..., nodes[end]..., weights_high..., weights_low... ]
function pack(
    nodes::Vector{SVector{1,T}}, weights_high::Vector{T}, weights_low::Vector{T}
) where {T}
    U = Vector{T}()
    for node in nodes
        append!(U, node)
    end
    append!(U, weights_high)
    append!(U, weights_low)
    return U
end

function unpack(U::Vector{T}, matrices::Vector{Vector{SMatrix{1,1,T,1}}}) where {T}
    H = length(matrices)

    nodes = [SVector{1,T}(t) for t in U[1:H]]
    weights_high = U[(H + 1):(2H)]
    weights_low = U[(2H + 1):end]

    return (
        matrices=matrices, nodes=nodes, weights_high=weights_high, weights_low=weights_low
    )
end

function increase_precision(quad)
    matrices = quad[:orbit_matrices]
    nodes = quad[:nodes]
    weights_high = quad[:weights_high]
    weights_low = quad[:weights_low]

    H, L = length(weights_high), length(weights_low)

    U = pack(nodes, weights_high, weights_low)

    return unpack(U, matrices)
end

const segment = hai.segment(big"-1.0", big"1.0")
const reference = hai.reference_domain(typeof(segment))

const matrix_0::Vector{AffineMap{1,BigFloat,1}} = [
    AffineMap(SMatrix{1,1}(big"1.0"), SVector{1}(big"0.0"))
]
const matrix_x::Vector{AffineMap{1,BigFloat,1}} = [
    AffineMap(SMatrix{1,1}(big"1.0"), SVector{1}(big"0.0")),
    AffineMap(SMatrix{1,1}(big"-1.0"), SVector{1}(big"0.0")),
]

# Values from QuadGK.jl
const gk7 = assemble(
    (matrix_0, "0.0000000000000000", "0.45091653865847414", "0.8888888888888885"),
    (matrix_x, "0.7745966692414834", "0.26848808986833345", "0.5555555555555556"),
    (matrix_x, "0.4342437493468026", "0.40139741477596220", ""),
    (matrix_x, "0.9604912687080203", "0.10465622602646725", ""),
)

const gk15 = assemble(
    (matrix_0, "0.0000000000000000", "0.20948214108472793", "0.41795918367346907"),
    (matrix_x, "0.4058451513773972", "0.19035057806478559", "0.38183005050511887"),
    (matrix_x, "0.7415311855993945", "0.14065325971552592", "0.27970539148927670"),
    (matrix_x, "0.9491079123427585", "0.06309209262997842", "0.12948496616886981"),
    (matrix_x, "0.2077849550078984", "0.20443294007529877", ""),
    (matrix_x, "0.5860872354676911", "0.16900472663926788", ""),
    (matrix_x, "0.8648644233597691", "0.10479001032225017", ""),
    (matrix_x, "0.9914553711208126", "0.02293532201052925", ""),
)

const gk31 = assemble(
    (matrix_0, "0.0000000000000000", "0.10133000701479164", "0.20257824192556112"),
    (matrix_x, "0.2011940939974345", "0.09917359872179185", "0.19843148532711144"),
    (matrix_x, "0.3941513470775634", "0.09312659817082533", "0.18616100001556220"),
    (matrix_x, "0.5709721726085388", "0.08308050282313299", "0.16626920581699406"),
    (matrix_x, "0.7244177313601701", "0.06985412131872824", "0.13957067792615430"),
    (matrix_x, "0.8482065834104272", "0.05348152469092807", "0.10715922046717188"),
    (matrix_x, "0.9372733924007060", "0.03534636079137582", "0.07036604748810820"),
    (matrix_x, "0.9879925180204855", "0.01500794732931611", "0.03075324199611720"),
    (matrix_x, "0.1011420669187175", "0.10076984552387552", ""),
    (matrix_x, "0.2991800071531688", "0.09664272698362363", ""),
    (matrix_x, "0.4850818636402397", "0.08856444305621179", ""),
    (matrix_x, "0.6509967412974169", "0.07684968075772045", ""),
    (matrix_x, "0.7904185014424660", "0.06200956780067060", ""),
    (matrix_x, "0.8972645323440819", "0.04458975132476482", ""),
    (matrix_x, "0.9677390756791391", "0.02546084732671520", ""),
    (matrix_x, "0.9980022986933971", "0.00537747987292333", ""),
)

# for (name, cbt, n) in [("G7", g7, 34), ("K15", k15, 34), ("G15", g15, 33), ("K31", k31, 33)]
#     println(">> $name <<")

#     local Φ = hai.map_to_reference(segment)
#     local j =
#         hai.abs_det_jac(hai.reference_orthotope(BigFloat, 1)) / hai.abs_det_jac(segment)

#     fmt_node = Format("[\"%.$(n)e\"],")
#     for x in cbt[:nodes]
#         local v = Φ(SVector{1}([x]))[1]
#         println(format(fmt_node, round(v; sigdigits=n + 1, base=10)))
#     end
#     println()

#     fmt_weight = Format("\"%.$(n)e\",")
#     for w in cbt[:weights]
#         println(format(fmt_weight, round(j * w; sigdigits=n + 1, base=10)))
#     end
#     println()
# end
