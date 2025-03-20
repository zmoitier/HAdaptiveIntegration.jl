using LinearAlgebra

function generate_monomials(k)
    monomials = Tuple{Int,Int}[]
    for m in 0:k
        for n in 0:(k - m)
            push!(monomials, (m, n))
        end
    end
    return monomials
end

function exact_integral(m, n)
    m_fact = factorial(big(m))
    n_fact = factorial(big(n))
    denom_fact = factorial(big(m + n + 2))
    return BigFloat(m_fact * n_fact / denom_fact)
end

function unpack(z)
    N = div(length(z), 3)
    x = [z[3 * i - 2] for i in 1:N]
    y = [z[3 * i - 1] for i in 1:N]
    w = [z[3 * i] for i in 1:N]
    return x, y, w
end

function compute_F(z, monomials, exact_integrals)
    x, y, w = unpack(z)
    F = similar(z, length(monomials))
    for (j, (m, n)) in enumerate(monomials)
        sum_q = zero(eltype(z))
        for i in eachindex(x, y, w)
            sum_q += w[i] * x[i]^m * y[i]^n
        end
        F[j] = sum_q - exact_integrals[j]
    end
    return F
end

function compute_J(z, monomials)
    x, y, w = unpack(z)
    N = length(x)
    M = length(monomials)
    J = zeros(eltype(z), M, 3N)
    for (j, (m, n)) in enumerate(monomials)
        for i in 1:N
            xi = x[i]
            yi = y[i]
            wi = w[i]
            # Compute derivatives
            dw = xi^m * yi^n
            dx = m == 0 ? zero(wi) : wi * m * xi^(m - 1) * yi^n
            dy = n == 0 ? zero(wi) : wi * n * xi^m * yi^(n - 1)
            # Assign to Jacobian
            J[j, 3i - 2] = dx
            J[j, 3i - 1] = dy
            J[j, 3i] = dw
        end
    end
    return J
end

function refine_quadrature(z0, monomials, exact_integrals; max_iter=10, tol=1e-80)
    z = copy(z0)
    for iter in 1:max_iter
        @show iter
        F = compute_F(z, monomials, exact_integrals)
        current_norm = norm(F, Inf)
        println("Iteration $iter, norm(F) = $current_norm")
        if current_norm < tol
            println("Converged after $iter iterations.")
            break
        end
        J = compute_J(z, monomials)
        @show eltype(J), eltype(F)
        delta = J \ (-F)
        z += delta
    end
    return z
end

# Example usage for k=1 (adjust as needed)
k = 5
N = div((k + 1) * (k + 2), 6)  # Ensure 3N equals number of monomials

x = [
    3.333333333333333e-01
    7.974269853530872e-01
    1.012865073234563e-01
    1.012865073234563e-01
    5.971587178976980e-02
    4.701420641051151e-01
    4.701420641051151e-01
]

y = [
    3.333333333333333e-01
    1.012865073234563e-01
    7.974269853530872e-01
    1.012865073234563e-01
    4.701420641051151e-01
    5.971587178976980e-02
    4.701420641051151e-01
]

w = [
    1.125000000000000e-01,
    6.296959027241359e-02,
    6.296959027241359e-02,
    6.296959027241359e-02,
    6.619707639425308e-02,
    6.619707639425308e-02,
    6.619707639425308e-02,
]

# Convert to BigFloat
z0 = Vector{BigFloat}(undef, 3 * N)
for i in 1:N
    z0[3i - 2] = BigFloat(x[i])
    z0[3i - 1] = BigFloat(y[i])
    z0[3i] = BigFloat(w[i])
end

monomials = generate_monomials(k)
exact_integrals = [exact_integral(m, n) for (m, n) in monomials]

@assert 3N == length(monomials) "System not square: 3N must equal number of monomials."

# Refine the quadrature
z = refine_quadrature(z0, monomials, exact_integrals; max_iter=10, tol=1e-100)

# # Unpack and print refined values
# x_refined, y_refined, w_refined = unpack(z)
# println("Refined nodes and weights:")
# for i in 1:N
#     println("Node $i: x = $(x_refined[i]), y = $(y_refined[i]), weight = $(w_refined[i])")
# end
