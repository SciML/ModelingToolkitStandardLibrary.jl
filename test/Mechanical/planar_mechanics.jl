# using Revise
# using Plots
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Mechanical.PlanarMechanics

tspan = (0.0, 3.0)
g = -9.807

@testset "Free body" begin
    # https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Examples/FreeBody.mo
    m = 2
    I = 1
    @named body = Body(; m, I)
    @named model = ODESystem(Equation[],
        t,
        [],
        [],
        systems = [body])
    sys = structural_simplify(model)
    unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan, []; jac = true)

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)

    free_falling_displacement = 0.5 * g * tspan[end]^2  # 0.5 * g * t^2
    @test sol[body.ry][end] ≈ free_falling_displacement
    @test sol[body.rx][end] == 0  # no horizontal displacement
    @test all(sol[body.phi] .== 0)
    # plot(sol, idxs = [body.rx, body.ry])
end

@testset "Pendulum" begin
    # https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Examples/Pendulum.mo
    @named ceiling = Fixed()
    @named rod = FixedTranslation(rx = 1.0, ry = 0.0)
    @named body = Body(m = 1, I = 0.1)
    @named revolute = Revolute()

    connections = [
        connect(ceiling.frame, revolute.frame_a),
        connect(revolute.frame_b, rod.frame_a),
        connect(rod.frame_b, body.frame)
    ]

    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [body, revolute, rod, ceiling])
    sys = structural_simplify(model)

    @test_broken length(unknowns(sys)) == 2
    unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan, []; jac = true)

    sol = solve(prob, Rodas5P())
    @test_broken SciMLBase.successful_retcode(sol)
end

@testset "Prismatic" begin
    # just testing instantiation
    @test_nowarn @named prismatic = Prismatic(x = 1.0, y = 0.0)
end

# @testset "AbsoluteAccCentrifugal" begin
#     # https://github.com/dzimmer/PlanarMechanics/blob/443b007bcc1522bb172f13012e2d7a8ecc3f7a9b/PlanarMechanicsTest/Sensors.mo#L221-L332
#     m = 1
#     I = 0.1
#     ω = 10
#     resolve_in_frame = :world

#     # components
#     @named body = Body(; m, I, gy = 0.0)
#     @named fixed_translation = FixedTranslation(; rx = 10.0, ry = 0.0)
#     @named fixed = Fixed()
#     @named revolute = Revolute(constant_ω = ω)

#     # sensors
#     @named abs_v_sensor = AbsoluteVelocity(; resolve_in_frame)

#     eqs = [
#         connect(fixed.frame, revolute.frame_a),
#         connect(revolute.frame_b, fixed_translation.frame_a),
#         connect(fixed_translation.frame_b, body.frame),
#         connect_sensor(body.frame, abs_v_sensor.frame_a)...
#     ]

#     @named model = ODESystem(eqs,
#         t,
#         [],
#         [],
#         systems = [
#             body,
#             fixed_translation,
#             fixed,
#             revolute,
#             abs_v_sensor
#         ])
#     sys = structural_simplify(model)
#     u0 = [0.0, ω, 0.0]
#     prob = ODEProblem(sys, u0, tspan, []; jac = true)
#     sol = solve(prob, Rodas5P())

#     # phi 
#     @test sol[body.phi][end] ≈ tspan[end] * ω
#     @test all(sol[body.ω] .≈ ω)

#     test_points = [i / ω for i in 0:0.1:10]

#     # instantaneous linear velocity
#     v_singal(t) = -ω^2 * sin.(ω .* t)
#     @test all(v_singal.(test_points) .≈ sol.(test_points; idxs = abs_v_sensor.v_x.u))

#     # instantaneous linear acceleration
#     a_singal(t) = -ω^3 * cos.(ω .* t)
#     @test all(a_singal.(test_points) .≈ sol.(test_points; idxs = body.ax))
# end

# @testset "Sensors (two free falling bodies)" begin
#     m = 1
#     I = 1
#     resolve_in_frame = :world

#     @named body1 = Body(; m, I)
#     @named body2 = Body(; m, I)
#     @named base = Fixed()

#     @named abs_pos_sensor = AbsolutePosition(; resolve_in_frame)
#     @named abs_v_sensor = AbsoluteVelocity(; resolve_in_frame)
#     @named abs_a_sensor = AbsoluteAcceleration(; resolve_in_frame)
#     @named rel_pos_sensor1 = RelativePosition(; resolve_in_frame)
#     @named rel_pos_sensor2 = RelativePosition(; resolve_in_frame)
#     @named rel_v_sensor1 = RelativeVelocity(; resolve_in_frame)
#     @named rel_v_sensor2 = RelativeVelocity(; resolve_in_frame)
#     @named rel_a_sensor1 = RelativeAcceleration(; resolve_in_frame)
#     @named rel_a_sensor2 = RelativeAcceleration(; resolve_in_frame)

#     connections = [
#         connect_sensor(body1.frame, abs_pos_sensor.frame_a)...,
#         connect_sensor(body1.frame, abs_v_sensor.frame_a)...,
#         connect_sensor(body1.frame, abs_a_sensor.frame_a)...,
#         connect_sensor(body1.frame, rel_pos_sensor1.frame_a)...,
#         connect_sensor(base.frame, rel_pos_sensor1.frame_b)...,
#         connect_sensor(body1.frame, rel_pos_sensor2.frame_a)...,
#         connect_sensor(body2.frame, rel_pos_sensor2.frame_b)...,
#         connect_sensor(base.frame, rel_v_sensor1.frame_a)...,
#         connect_sensor(body1.frame, rel_v_sensor1.frame_b)...,
#         connect_sensor(body1.frame, rel_v_sensor2.frame_a)...,
#         connect_sensor(body2.frame, rel_v_sensor2.frame_b)...,
#         connect_sensor(body1.frame, rel_a_sensor1.frame_a)...,
#         connect_sensor(base.frame, rel_a_sensor1.frame_b)...,
#         connect_sensor(body1.frame, rel_a_sensor2.frame_a)...,
#         connect_sensor(body2.frame, rel_a_sensor2.frame_b)...
#     ]

#     @named model = ODESystem(connections,
#         t,
#         [],
#         [],
#         systems = [
#             body1,
#             body2,
#             base,
#             abs_pos_sensor,
#             abs_v_sensor,
#             abs_a_sensor,
#             rel_pos_sensor1,
#             rel_pos_sensor2,
#             rel_v_sensor1,
#             rel_v_sensor2,
#             rel_a_sensor1,
#             rel_a_sensor2
#         ])

#     sys = structural_simplify(model)
#     unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
#     prob = ODEProblem(sys, unset_vars .=> 0.0, tspan, []; jac = true)

#     sol = solve(prob, Rodas5P())
#     @test SciMLBase.successful_retcode(sol)

#     # the two bodyies falled the same distance, and so the absolute sensor attached to body1
#     @test sol[abs_pos_sensor.y.u][end] ≈ sol[body1.ry][end] ≈ sol[body2.ry][end] ≈
#           0.5 * g * tspan[end]^2

#     # sensor1 is attached to body1, so the relative y-position between body1 and the base is
#     # equal to the absolute y-position of body1
#     @test sol[body1.ry][end] ≈ -sol[rel_pos_sensor1.rel_y.u][end]

#     # the relative y-position between body1 and body2 is zero
#     @test sol[rel_pos_sensor2.rel_y.u][end] == 0

#     # no displacement in the x-direction
#     @test sol[abs_pos_sensor.x.u][end] ≈ sol[body1.rx][end] ≈ sol[body2.rx][end]

#     # velocity after t seconds v = g * t, so the relative y-velocity between body1 and the base is
#     # equal to the absolute y-velocity of body1
#     @test sol[abs_v_sensor.v_y.u][end] ≈ sol[rel_v_sensor1.rel_v_y.u][end] ≈ g * tspan[end]

#     # the relative y-velocity between body1 and body2 is zero
#     @test sol[rel_v_sensor2.rel_v_y.u][end] == 0

#     # the body is under constant acclertation = g
#     @test all(sol[abs_a_sensor.a_y.u] .≈ g)

#     # the relative y-accleration between body1 and the base is
#     # equal to the absolute y-accleration of body1
#     @test sol[abs_a_sensor.a_y.u][end] ≈ -sol[rel_a_sensor1.rel_a_y.u][end]

#     # the relative y-accleration between body1 and body2 is zero
#     @test sol[rel_a_sensor2.rel_a_y.u][end] == 0
# end

# @testset "Measure Demo" begin
#     # https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Examples/MeasureDemo.mo
#     @named body = Body(; m = 1, I = 0.1)
#     @named fixed_translation = FixedTranslation(; rx = 1, ry = 0)
#     @named fixed = Fixed()
#     @named body1 = Body(; m = 0.4, I = 0.02)
#     @named fixed_translation1 = FixedTranslation(; rx = 0.4, ry = 0)
#     @named abs_pos_sensor = AbsolutePosition(; resolve_in_frame = :world)
#     @named rel_pos_sensor = RelativePosition(; resolve_in_frame = :world)
#     @named revolute1 = Revolute()
#     @named abs_v_sensor = AbsoluteVelocity(; resolve_in_frame = :frame_a)
#     @named rel_v_sensor = RelativeVelocity(; resolve_in_frame = :frame_b)
#     @named abs_a_sensor = AbsoluteAcceleration(; resolve_in_frame = :world)
#     @named rel_a_sensor = RelativeAcceleration(; resolve_in_frame = :frame_b)
#     @named revolute2 = Revolute()

#     connections = [
#         connect(fixed_translation.frame_b, body.frame),
#         connect(fixed_translation1.frame_b, body1.frame),
#         connect(fixed.frame, revolute1.frame_a),
#         connect(revolute1.frame_b, fixed_translation.frame_a),
#         # connect(abs_a_sensor.frame_resolve, abs_a_sensor.frame_a),
#         connect(revolute2.frame_b, fixed_translation1.frame_a),
#         connect(revolute2.frame_a, fixed_translation.frame_b),
#         # connect_sensor(fixed_translation.frame_b, rel_a_sensor.frame_a)...,
#         # connect(fixed_translation.frame_b, rel_v_sensor.frame_a),
#         # connect(fixed_translation.frame_b, rel_v_sensor.frame_a),
#         # connect(rel_a_sensor.frame_b, body1.frame_a),
#         # connect(rel_v_sensor.frame_b, body1.frame_a),
#         # connect(rel_v_sensor.frame_b, body1.frame_a),
#         connect_sensor(body1.frame, abs_a_sensor.frame_a)...        # connect_sensor(body1.frame, abs_v_sensor.frame_a)...,        # connect_sensor(body1.frame, abs_pos_sensor.frame_a)...,
#     ]

#     @named model = ODESystem(connections,
#         t,
#         [],
#         [],
#         systems = [
#             fixed_translation,
#             body,
#             fixed,
#             body1,
#             fixed_translation1,
#             revolute1,
#             revolute2,
#             abs_pos_sensor
#         ])
#     sys = structural_simplify(model)
#     unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
#     prob = ODEProblem(sys, unset_vars .=> 0.0, (0, 5), []; jac = true)
#     sol = solve(prob, Rodas5P())
#     @test_broken SciMLBase.successful_retcode(sol)
# end

# @testset "SpringDamper" begin
#     # https://github.com/dzimmer/PlanarMechanics/blob/master/PlanarMechanics/Examples/SpringDamperDemo.mo
#     @named spring_damper = SpringDamper(;
#         s_relx0 = 0,
#         d_y = 1,
#         s_rely0 = 0,
#         d_phi = 0,
#         c_y = 5,
#         c_x = 5,
#         d_x = 1,
#         c_phi = 0)
#     @named body = Body(; I = 0.1, m = 0.5, rx = 1, ry = 1)
#     @named fixed = Fixed()
#     @named fixed_translation = FixedTranslation(; rx = -1, ry = 0)

#     connections = [
#         connect(fixed.frame, fixed_translation.frame_a),
#         connect(fixed_translation.frame_b, spring_damper.frame_a),
#         connect(spring_damper.frame_b, body.frame)
#     ]
#     @named model = ODESystem(connections,
#         t,
#         [],
#         [],
#         systems = [
#             spring_damper,
#             body,
#             fixed,
#             fixed_translation
#         ])
#     sys = structural_simplify(model)
#     unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
#     prob = ODEProblem(sys, unset_vars .=> 0.0, (0, 5), []; jac = true)
#     sol = solve(prob, Rodas5P())
#     @test SciMLBase.successful_retcode(sol)
# end

# @testset "Spring and damper demo" begin
#     # https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Examples/SpringDemo.mo
#     @named body = Body(; m = 0.5, I = 0.1)
#     @named fixed = Fixed()
#     @named spring = Spring(; c_y = 10, s_rely0 = -0.5, c_x = 1, c_phi = 1e5)
#     @named damper = Damper(d = 1)
#     @named prismatic = Prismatic(; x = 0, y = 1)

#     connections = [
#         connect(fixed.frame, spring.frame_a),
#         connect(spring.frame_b, body.frame),
#         connect(damper.frame_a, spring.frame_a),
#         connect(damper.frame_b, spring.frame_b),
#         connect(spring.frame_a, prismatic.frame_a),
#         connect(prismatic.frame_b, spring.frame_b)
#     ]

#     @named model = ODESystem(connections,
#         t,
#         [],
#         [],
#         systems = [
#             body,
#             fixed,
#             spring,
#             damper,
#             prismatic
#         ])
#     sys = structural_simplify(model)
#     unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
#     prob = ODEProblem(sys, unset_vars .=> 0.0, (0, 5), []; jac = true)
#     sol = solve(prob, Rodas5P())
#     @test SciMLBase.successful_retcode(sol)
# end
