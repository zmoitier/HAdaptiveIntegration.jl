import HAdaptiveIntegration as hai

using Printf
using StaticArrays

D = 3
gm = hai.GrundmannMoeller{D}(9)
ec = hai.embedded_cubature(Rational{Int}, gm)
L = length(ec.weights_low)

# println("-"^44)
# for (node, wh, wl) in zip(ec.nodes, ec.weights_high, ec.weights_low)
#     x, y = node
#     println(@sprintf("(%.4e, %.4e) %+.2e %+.2e", x, y, wh, wl))
# end
# for (node, wh) in zip(ec.nodes[(L + 1):end], ec.weights_high[(L + 1):end])
#     x, y = node
#     println(@sprintf("(%.4e, %.4e) %+.2e", x, y, wh))
# end

counter = Dict{SVector{D,Rational{Int}},Int}()
for node in ec.nodes
    counter[node] = get(counter, node, 0) + 1
end
for (node, count) in counter
    if count > 1
        println(@sprintf("%s: %d", node, count))
    end
end

display(sum(ec.weights_high))
display(sum(ec.weights_low))
