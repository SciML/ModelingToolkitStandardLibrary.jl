using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Mechanical.PlanarMechanics

@parameters t
D = Differential(t)

@testset "Fixed" begin
    @test true
end
