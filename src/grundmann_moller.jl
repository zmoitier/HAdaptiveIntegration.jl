function grundmann_moller(n::Int, s::Int)
    s_prime = s - 1  # Grundmann-Möller parameter for degree 2s+1 = 2(s-1)+1 = 2s-1
    m = 2 * s_prime + n + 1  # Common denominator for barycentric coordinates
    points = Vector{Vector{Float64}}()
    weights = Vector{Float64}()

    for k in 0:s_prime
        # Skip if t_sum = s_prime - k is negative
        t_sum = s_prime - k
        t_sum < 0 && continue

        # Generate all combinations of (k+1) coordinates from (n+1) vertices
        for subset in combinations(n + 1, k + 1)
            # Generate all non-negative integer tuples (t_0, ..., t_k) summing to t_sum
            ts = multiexponents(k + 1, t_sum)
            for t in ts
                bary = fill(1.0 / m, n + 1)  # Default value for unselected coordinates
                # Update selected coordinates with (2t_i + 1)/m
                for (i, idx) in enumerate(subset)
                    bary[idx] = (2 * t[i] + 1) / m
                end
                # Check barycentric sum
                @assert isapprox(sum(bary), 1.0, atol=1e-10) "Barycentric sum ≠ 1"
                push!(points, bary)
                # Compute weight (approximate formula; may need adjustment)
                weight_num = (-1.0)^(s_prime - k) * factorial(k)^2
                weight_den = factorial(s_prime + n + 1) * prod(2 * t_i + 1 for t_i in t)
                weight = (2.0^(-2 * s_prime)) * weight_num / weight_den
                push!(weights, weight)
            end
        end
    end

    # Convert barycentric to Cartesian coordinates (drop the first coordinate)
    cartesian = [p[2:end] for p in points]

    # Normalize weights to integrate 1 over the simplex (volume = 1/n!)
    volume = 1.0 / factorial(n)
    total_weight = sum(weights)
    weights .= (weights ./ total_weight) .* volume

    return (cartesian, weights)
end

function multiexponents(dim::Int, sum_total::Int)
    exponents = Vector{Int}[]
    if dim == 1
        return [[sum_total]]
    else
        for i in 0:sum_total
            for rest in multiexponents(dim - 1, sum_total - i)
                push!(exponents, [i; rest])
            end
        end
    end
    return exponents
end
