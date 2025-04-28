using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using Quadmath
using Test

function validate_orders(
    ec::EmbeddedCubature{D,T},
    get_exact_values::Function,
    order_high::Int,
    order_low::Int;
    rtol=10 * eps(T),
) where {D,T}
    exact_values = get_exact_values(D, order_high)

    for k in 0:order_high
        for (α, vₑₓ) in exact_values[k + 1]
            fct = x -> prod(xᵢ^eᵢ for (xᵢ, eᵢ) in zip(x, α))
            Iₕ, Iₗ = ec(fct)

            if k ≤ order_low
                if !isapprox(Iₗ, vₑₓ; rtol=rtol)
                    @error "fail to integrate the low order cubature within tolerance at degree = $α. Relative error: $(abs(Iₗ/vₑₓ-1)), rtol: $rtol."
                    return false
                end
            end

            if !isapprox(Iₕ, vₑₓ; rtol=rtol)
                @error "fail to integrate the high order cubature within tolerance at degree = $α. Relative error: $(abs(Iₕ/vₑₓ-1)), rtol: $rtol."
                return false
            end
        end
    end

    return true
end

function (ec::EmbeddedCubature{D,T})(fct) where {D,T}
    H, L = length(ec.weights_high), length(ec.weights_low)

    v = fct(ec.nodes[1])
    Iₗₒ = ec.weights_low[1] * v
    Iₕᵢ = ec.weights_high[1] * v
    for i in 2:L
        v = fct(ec.nodes[i])
        Iₗₒ += ec.weights_low[i] * v
        Iₕᵢ += ec.weights_high[i] * v
    end
    for i in (L + 1):H
        Iₕᵢ += ec.weights_high[i] * fct(ec.nodes[i])
    end
    return Iₕᵢ, Iₗₒ
end

function integral_monomial_simplex(dim::Int, k_max::Int)
    @assert (dim > 0) && (k_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    mlt_idx_val = [[(i,) => 1//prod((i + 1):(i + dim))] for i in 0:k_max]

    for d in 2:dim
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:k_max]
        for (k, mi_val) in zip(Iterators.countfrom(0), mlt_idx_val)
            for (mi, val) in mi_val
                push!(new[k + 1], (0, mi...) => val)
                v = 1
                for n in 1:(k_max - k)
                    v *= n//(n + k + dim)
                    push!(new[k + 1 + n], (n, mi...) => val * v)
                end
            end
        end
        mlt_idx_val = new
    end

    return mlt_idx_val
end

function integral_monomial_orthotope(dim::Int, k_max::Int)
    @assert (dim > 0) && (k_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    mlt_idx_val = [[(i,) => 1//(i + 1)] for i in 0:k_max]

    for d in 2:dim
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:k_max]
        for (k, mi_val) in zip(Iterators.countfrom(0), mlt_idx_val)
            for (mi, val) in mi_val
                for n in 0:(k_max - k)
                    push!(new[k + 1 + n], (n, mi...) => val//(n + 1))
                end
            end
        end
        mlt_idx_val = new
    end

    return mlt_idx_val
end

@testset "Rule construction" begin
    @testset "Tabulated embedded cubature" begin
        tec = TabulatedEmbeddedCubature{Orthotope{1}}(;
            description="Gauss (SEGMENT_G3)",
            reference="",
            precision=16,
            nodes=[["5e-1"], ["1.127016653792583e-1"], ["8.872983346207417e-1"]],
            weights_high=[
                "4.444444444444444e-1", "2.777777777777778e-1", "2.777777777777778e-1"
            ],
            order_high=5,
            weights_low=["1"],
            order_low=1,
        )
        @test typeof(tec) <: TabulatedEmbeddedCubature{Orthotope{1}}
        @test orders(tec) == (5, 1)

        ec = embedded_cubature(tec)
        @test typeof(ec) <: EmbeddedCubature{1,Float64}

        ec_ref = embedded_cubature(
            [[0.5], [(1 - √(3 / 5)) / 2], [(1 + √(3 / 5)) / 2]],
            [4 / 9, 5 / 18, 5 / 18],
            [1.0],
        )
        @test typeof(ec_ref) <: EmbeddedCubature{1,Float64}

        @test ec.nodes ≈ ec_ref.nodes
        @test ec.weights_high ≈ ec_ref.weights_high
        @test ec.weights_low ≈ ec_ref.weights_low
    end

    @testset "Radon-Laurie" begin
        rl = RadonLaurie()
        @test typeof(rl) <: AbstractRule{Simplex{2}}

        @test orders(rl) == (8, 5)
        @test validate_orders(embedded_cubature(rl), integral_monomial_simplex, 8, 5)
    end

    @testset "Grundmann-Möller" begin
        @test_throws AssertionError typeof(GrundmannMoeller{1}(7, 6))
        @test_throws AssertionError typeof(GrundmannMoeller{1}(5, 7))

        gm = GrundmannMoeller{2}(5, 3)
        @test typeof(gm) <: AbstractRule{Simplex{2}}
        @test orders(gm) == (5, 3)
        @test validate_orders(embedded_cubature(gm), integral_monomial_simplex, 5, 3)
    end

    @testset "Genz-Malik" begin
        gm = GenzMalik{2}()
        @test typeof(gm) <: AbstractRule{Orthotope{2}}
        @test orders(gm) == (7, 5)
        @test validate_orders(embedded_cubature(gm), integral_monomial_orthotope, 7, 5)
    end
end

@testset "Tabulated rules" begin
    T = Float128

    @testset "Segment" begin
        ec = embedded_cubature(GrundmannMoeller{1}(7, 5), T)
        @test validate_orders(ec, integral_monomial_simplex, 7, 5)

        ec = embedded_cubature(GenzMalik{1}(), T)
        @test validate_orders(ec, integral_monomial_orthotope, 7, 5)
    end

    @testset "Triangle" begin
        for tec in (TRIANGLE_GM19, TRIANGLE_RL19)
            ec = embedded_cubature(tec, T)
            @test validate_orders(ec, integral_monomial_simplex, orders(tec)...)
        end
    end

    @testset "Square" begin
        for tec in (SQUARE_GM17, SQUARE_CH21, SQUARE_CH25)
            ec = embedded_cubature(tec, T)
            @test validate_orders(ec, integral_monomial_orthotope, orders(tec)...)
        end
    end

    @testset "Tetrahedron" begin
        tec = TETRAHEDRON_GM35
        ec = embedded_cubature(tec, T)
        @test validate_orders(ec, integral_monomial_simplex, orders(tec)...)
    end

    @testset "Cube" begin
        for tec in (CUBE_GM33, CUBE_BE65, CUBE_BE115)
            ec = embedded_cubature(tec, T)
            @test validate_orders(ec, integral_monomial_orthotope, orders(tec)...)
        end
    end
end
