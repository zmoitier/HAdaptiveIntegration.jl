using Test
import HAdaptiveIntegration as HAI

@testset "domain_subdivision.jl" begin
    @testset "segment" begin
        s = HAI.segment(-1.0, 2.0)

        s2 = HAI.subdivide_segment2(s)
        @test length(s2) == 2
        @test sum(HAI.abs_jacobian_determinant.(s2)) ≈ HAI.abs_jacobian_determinant(s)
    end

    @testset "triangle" begin
        s = HAI.triangle((-1.0, -2.0), (1.0, 0.0), (0.0, 2.0))

        s2 = HAI.subdivide_triangle2(s)
        @test length(s2) == 2
        @test sum(HAI.abs_jacobian_determinant.(s2)) ≈ HAI.abs_jacobian_determinant(s)

        s4 = HAI.subdivide_triangle4(s)
        @test length(s4) == 4
        @test sum(HAI.abs_jacobian_determinant.(s4)) ≈ HAI.abs_jacobian_determinant(s)
    end

    @testset "rectangle" begin
        s = HAI.rectangle((-2.0, 1.0), (1.0, 2.0))

        s2 = HAI.subdivide_rectangle4(s)
        @test length(s2) == 4
        @test sum(HAI.abs_jacobian_determinant.(s2)) ≈ HAI.abs_jacobian_determinant(s)
    end

    @testset "tetrahedron" begin
        s = HAI.tetrahedron(
            (-1.0, -2.0, -3.0), (1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)
        )

        s2 = HAI.subdivide_tetrahedron8(s)
        @test length(s2) == 8
        @test sum(HAI.abs_jacobian_determinant.(s2)) ≈ HAI.abs_jacobian_determinant(s)
    end

    @testset "cuboid" begin
        s = HAI.cuboid((-3.0, 1.0), (-1.0, 2.0), (-2.0, 3.0))

        s2 = HAI.subdivide_cuboid8(s)
        @test length(s2) == 8
        @test sum(HAI.abs_jacobian_determinant.(s2)) ≈ HAI.abs_jacobian_determinant(s)
    end
end
