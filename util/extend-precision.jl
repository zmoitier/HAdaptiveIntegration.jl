import HAdaptiveIntegration as HAI

function exponents(partial_degree_max::Int)
    αs = Vector{NTuple{2,Int}}()

    if partial_degree_max ≤ 0
        return αs
    end

    for α in Base.product(0:partial_degree_max, 0:partial_degree_max)
        push!(αs, α)
    end

    sort!(αs; by=α -> (sum(α), α[2], α[1]))
    return αs
end

setprecision(35; base=10)

square = HAI.rectangle(BigFloat.(["0.0", "0.0"]), BigFloat.(["1.0", "1.0"]))
ec = HAI.embedded_cubature(HAI.SQUARE_CH21G25, BigFloat)

nodes = ec.nodes
w_hi = ec.weights_high
w_lo = ec.weights_low

display(exponents(9))
