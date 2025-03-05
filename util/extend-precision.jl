import HAdaptiveIntegration as HAI

function exponents_even(partial_degree_max::Int)
    αs = Vector{NTuple{2,Int}}()

    if partial_degree_max ≤ 0
        return αs
    end

    for α in Base.product(0:partial_degree_max, 0:partial_degree_max)
        if (α[1] % 2 == 0) && (α[2] % 2 == 0)
            push!(αs, α)
        end
    end

    sort!(αs; by=α -> (sum(α), α[2], α[1]))
    return αs
end

setprecision(35; base=10)

square = HAI.rectangle(BigFloat.(["-1.0", "-1.0"]), BigFloat.(["1.0", "1.0"]))
ec = HAI.embedded_cubature(HAI.SQUARE_CH21G25, BigFloat)

Φ = HAI.map_from_reference(square)
J =
    HAI.abs_det_jacobian(square) /
    HAI.abs_det_jacobian(HAI.reference_domain(HAI.Rectangle{BigFloat}))

nodes = Φ.(ec.nodes)
w_hi = ec.weights_high .* J
w_lo = ec.weights_low .* J

# display(exponents_even(8))
