using HAdaptiveIntegration.Domain:
    AbstractDomain,
    Cuboid,
    Orthotope,
    Rectangle,
    Segment,
    Simplex,
    Tetrahedron,
    Triangle,
    reference_domain
using HAdaptiveIntegration: integrate
using Test
using Unitful

# test for https://github.com/zmoitier/HAdaptiveIntegration.jl/issues/56

dimension(::Type{<:AbstractDomain{D}}) where {D} = D

@testset "Support for Unitful.jl" begin
    U = typeof(0.0u"m")
    T = typeof(one(U))

    @testset "dimensionless domain with unit integrand" begin
        for dom_type in (
            Segment{T},
            Triangle{T},
            Rectangle{T},
            Tetrahedron{T},
            Cuboid{T},
            Simplex{4,T},
            Orthotope{4,T},
        )
            domain = reference_domain(dom_type)
            I, E = @inferred integrate(p -> sum(p) * 1u"g", domain)
            @test unit(I) == u"g"
            @test unit(E) == u"g"
        end
    end

    @testset "unit domain with dimensionless integrand" begin
        for dom_type in (
            Segment{U},
            Triangle{U},
            Rectangle{U},
            Tetrahedron{U},
            Cuboid{U},
            Simplex{4,U},
            Orthotope{4,U},
        )
            domain = reference_domain(dom_type)
            dim = dimension(dom_type)
            I, E = @inferred integrate(p -> 1, domain)
            @test unit(I) == unit(1u"m"^dim)
            @test unit(E) == unit(1u"m"^dim)
        end
    end

    @testset "unit domain with unit integrand" begin
        for dom_type in (
            Segment{U},
            Triangle{U},
            Rectangle{U},
            Tetrahedron{U},
            Cuboid{U},
            Simplex{4,U},
            Orthotope{4,U},
        )
            domain = reference_domain(dom_type)
            dim = dimension(dom_type)
            I, E = @inferred integrate(p -> 1u"g", domain)
            @test unit(I) == unit(1u"g" * 1u"m"^dim)
            @test unit(E) == unit(1u"g" * 1u"m"^dim)
        end
    end
end
