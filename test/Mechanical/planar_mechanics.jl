using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Mechanical.PlanarMechanics
# using Plots

@parameters t
D = Differential(t)
tspan = (0.0, 3.0)

@testset "Free body" begin
    @named body = Body(m = 1, j = 0.1)
    @named model = ODESystem(Equation[],
        t,
        [],
        [],
        systems = [body])
    sys = structural_simplify(model)
    unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan, []; jac = true)

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)

    free_falling_displacement = -0.5 * 9.807 * tspan[2]^2  # 0.5 * g * t^2
    @test sol[body.ry][end] ≈ free_falling_displacement
    @test sol[body.rx][end] == 0  # no horizontal displacement
    # plot(sol, idxs = [body.rx, body.ry])
end

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
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan, []; jac = true)
    sol = solve(prob, Rodas5P())

    # phi and omega for the pendulum body
    @test length(states(sys)) == 2
end

@testset "Prismatic" begin
    r = [1.0, 0.0]
    e = r / sqrt(r' * r)
    @named prismatic = Prismatic(rx = r[1], ry = r[2], ex = e[1], ey = e[2])
    # just testing instantiation
    @test true
end
