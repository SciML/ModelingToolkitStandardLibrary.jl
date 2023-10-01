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

@testset "Position Sensors (two free falling bodies)" begin
    m = 1
    j = 0
    resolve_in_frame = :world

    @named body1 = Body(; m, j)
    @named body2 = Body(; m, j)
    @named base = Fixed()

    @named abs_pos_sensor = AbsolutePosition(; resolve_in_frame)
    @named rel_pos_sensor1 = RelativePosition(; resolve_in_frame)
    @named rel_pos_sensor2 = RelativePosition(; resolve_in_frame)

    connections = [
        connect(body1.frame, abs_pos_sensor.frame_a),
        connect(rel_pos_sensor1.frame_a, body1.frame),
        connect(rel_pos_sensor1.frame_b, base.frame),
        connect(rel_pos_sensor2.frame_b, body1.frame),
        connect(rel_pos_sensor2.frame_a, body2.frame),
        body1.phi ~ 0,
        body1.fx ~ 0,
        body1.fy ~ m * -9.807,
        body2.phi ~ 0,
        body2.fx ~ 0,
        body2.fy ~ m * -9.807,
    ]

    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [body1, body2, base, abs_pos_sensor, rel_pos_sensor1, rel_pos_sensor2])

    sys = structural_simplify(model)
    unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan, []; jac = true)

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)

    # the two bodyies falled the same distance, and so the absolute sensor attached to body1
    @test sol[abs_pos_sensor.y.u][end] ≈ sol[body1.ry][end] ≈ sol[body2.ry][end]

    # sensor1 is attached to body1, so the relative y-position between body1 and the base is
    # equal to the y-position of body1
    @test sol[body1.ry][end] ≈ -sol[rel_pos_sensor1.rel_y.u][end]

    # the relative y-position between body1 and body2 is zero
    @test sol[rel_pos_sensor2.rel_y.u][end] == 0

    # no displacement in the x-direction
    @test sol[abs_pos_sensor.x.u][end] ≈ sol[body1.rx][end] ≈ sol[body2.rx][end]
end

@testset "Measure Demo" begin
    @test_nowarn @named abs_pos_w = AbsolutePosition(; resolve_in_frame = :world)
    @test_nowarn @named abs_pos_fa = AbsolutePosition(; resolve_in_frame = :frame_a)
    @test_nowarn @named abs_pos_fr = AbsolutePosition(; resolve_in_frame = :frame_resolve)

    @test_nowarn @named rel_pos_w = RelativePosition(; resolve_in_frame = :world)
    @test_nowarn @named rel_pos_fa = RelativePosition(; resolve_in_frame = :frame_a)
    @test_nowarn @named rel_pos_fb = RelativePosition(; resolve_in_frame = :frame_b)
    @test_nowarn @named rel_pos_fr = RelativePosition(; resolve_in_frame = :frame_resolve)
end
