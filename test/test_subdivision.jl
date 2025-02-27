using Test
import HAdaptiveIntegration as HAI

@testset "domain_subdivision.jl" begin
    @testset "segment" begin
        s = HAI.segment(-0.81, 0.48)
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_segment2))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_segment3))
    end

    @testset "triangle" begin
        s = HAI.triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_triangle2))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_triangle4))
    end

    @testset "rectangle" begin
        s = HAI.rectangle((-0.07, 0.42), (0.35, 0.71))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_rectangle4))
    end

    @testset "tetrahedron" begin
        s = HAI.tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_tetrahedron8))
    end

    @testset "cuboid" begin
        s = HAI.cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_cuboid8))
    end
end
