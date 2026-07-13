using ModelingToolkitBase: System, connect, mtkcompile, @named, t_nounits
using ModelingToolkitStandardLibrary.Blocks: Constant, StateSpace
using SciMLBase: ODEProblem
using Test

@testset "Array state codegen" begin
    A = [0 1; -1 -0.5]
    B = [0, 1]
    C = [0.9 1;]
    D = [0;;]
    @named ss = StateSpace(; A, B, C, D, x = zeros(2))
    @named source = Constant(; k = 1)
    @named model = System(
        connect(source.output, ss.input),
        t_nounits,
        systems = [ss, source]
    )

    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    @test prob.u0 == zeros(2)
end
