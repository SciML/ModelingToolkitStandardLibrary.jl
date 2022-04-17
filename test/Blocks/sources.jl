using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks

@parameters t

@testset "Sources" begin
    sources = [Constant, Sine, Cosine, Clock, Ramp, Step, ExpSine]
    for source in sources
    @named src = source()
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    @test_nowarn sol = solve(prob, Rodas4())
end