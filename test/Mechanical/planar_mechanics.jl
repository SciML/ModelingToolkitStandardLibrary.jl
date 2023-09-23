using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Mechanical.PlanarMechanics

@parameters t
D = Differential(t)

@testset "Pendulum" begin
    @named ceiling = Fixed()
    @named rod = FixedTranslation(rx = 1.0, ry = 0.0)
    @named body = Body(m = 1, j = 0.1)
    @named revolute = Revolute(phi = 0.0, ω = 0.0)

    connections = [
        connect(ceiling.frame, revolute.frame_a),
        connect(revolute.frame_b, rod.frame_a),
        connect(rod.frame_b, body.frame),
    ]

    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [body, revolute, rod, ceiling])
    sys = structural_simplify(model)
    unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, (0.0, 60), []; jac = true)
    sol = solve(prob, Rodas5P())

    # phi and omega for the pendulum body
    @test length(states(sys)) == 2
end
