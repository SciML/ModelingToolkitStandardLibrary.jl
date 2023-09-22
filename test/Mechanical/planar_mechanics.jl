using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Mechanical.PlanarMechanics

@parameters t
D = Differential(t)

@testset "Pendulum" begin
    @named body = Body(m = 1, j = 0.1)
    @named revolute = Revolute(phi = 0.0, ω = 0.0)
    @named rod = FixedTranslation(rx = 1.0, ry = 0.0)
    @named ground = Fixed()
    @test true
end
