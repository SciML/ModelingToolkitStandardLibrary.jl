using ModelingToolkitStandardLibrary.Mechanical.Rotational,
      ModelingToolkit, OrdinaryDiffEq,
      Test
using ModelingToolkit: t, D
import ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq: ReturnCode.Success

# using Plots

@testset "two inertias" begin
    @named fixed = Fixed()
    @named inertia1 = Inertia(J = 2) # this one is fixed
    @named spring = Spring(c = 1e4)
    @named damper = Damper(d = 10)
    @named inertia2 = Inertia(J = 2, phi = pi / 2)

    connections = [connect(fixed.flange, inertia1.flange_b)
                   connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
                   connect(spring.flange_b, damper.flange_b, inertia2.flange_a)]

    @named model = ODESystem(connections, t,
        systems = [fixed, inertia1, inertia2, spring, damper])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[], (0, 10.0))
    sol1 = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol1)

    prob = ODAEProblem(sys, Pair[], (0, 10.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)

    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, Pair[], (0, 10.0))
    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[inertia1.w] .== 0)
    @test sol[inertia2.w][end]≈0 atol=1e-3 # all energy has dissipated

    @named springdamper = SpringDamper(; c = 1e4, d = 10)
    connections = [connect(fixed.flange, inertia1.flange_b)
                   connect(inertia1.flange_b, springdamper.flange_a)
                   connect(springdamper.flange_b, inertia2.flange_a)]

    @named model = ODESystem(connections, t,
        systems = [fixed, inertia1, inertia2, springdamper])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[], (0, 10.0))
    sol2 = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)

    @test sol2(0:1:10, idxs = inertia2.w).u≈sol1(0:1:10, idxs = inertia2.w).u atol=1e-3

    # Plots.plot(sol; vars=[inertia1.w, inertia2.w])
end

@testset "two inertias with driving torque" begin
    amplitude = 10 # Amplitude of driving torque
    frequency = 5 # Frequency of driving torque
    J_motor = 0.1 # Motor inertia

    @named fixed = Fixed()
    @named torque = Torque(; use_support = true)
    @named inertia1 = Inertia(J = 2, phi = pi / 2)
    @named spring = Rotational.Spring(c = 1e4)
    @named damper = Damper(d = 10)
    @named inertia2 = Inertia(J = 4)
    @named sine = Blocks.Sine(amplitude = amplitude, frequency = frequency)

    connections = [connect(sine.output, torque.tau)
                   connect(torque.support, fixed.flange)
                   connect(torque.flange, inertia1.flange_a)
                   connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
                   connect(spring.flange_b, damper.flange_b, inertia2.flange_a)]

    @named model = ODESystem(connections, t,
        systems = [
            fixed,
            torque,
            inertia1,
            inertia2,
            spring,
            damper,
            sine
        ])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0,
        [D(D(inertia2.phi)) => 1.0; D.(states(model)) .=> 0.0], (0, 10.0))
    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)

    prob = ODAEProblem(sys, Pair[], (0, 1.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)

    @test all(isapprox.(sol[inertia1.w], -sol[inertia2.w] * 2, atol = 1)) # exact opposite oscillation with smaller amplitude J2 = 2*J1
    @test all(sol[torque.flange.tau] .== -sol[sine.output.u]) # torque source is equal to negative sine

    ## Test with constant torque source
    @named torque = ConstantTorque(use_support = true, tau_constant = 1)
    connections = [connect(torque.support, fixed.flange)
                   connect(torque.flange, inertia1.flange_a)
                   connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
                   connect(spring.flange_b, damper.flange_b, inertia2.flange_a)]

    @named model = ODESystem(connections, t,
        systems = [
            fixed,
            torque,
            inertia1,
            inertia2,
            spring,
            damper
        ])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[], (0, 10.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    @test sol(sol.t[end], idxs = inertia1.w)≈sol(sol.t[end], idxs = inertia2.w) rtol=0.1 # both inertias have same angular velocity after initial transient
end

# see: https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Mechanics.Rotational.Examples.First.html
@testset "first example" begin
    amplitude = 10 # Amplitude of driving torque
    frequency = 5 # Frequency of driving torque
    J_motor = 0.1 # Motor inertia
    J_load = 2 # Load inertia
    ratio = 10 # Gear ratio
    damping = 10 # Damping in bearing of gear

    @named fixed = Fixed()
    @named torque = Torque(use_support = true)
    @named inertia1 = Inertia(J = J_motor)
    @named idealGear = IdealGear(ratio = ratio, use_support = true)
    @named inertia2 = Inertia(J = 2)
    @named spring = Spring(c = 1e4)
    @named inertia3 = Inertia(J = J_load)
    @named damper = Damper(d = damping)
    @named sine = Blocks.Sine(amplitude = amplitude, frequency = frequency)

    connections = [connect(inertia1.flange_b, idealGear.flange_a)
                   connect(idealGear.flange_b, inertia2.flange_a)
                   connect(inertia2.flange_b, spring.flange_a)
                   connect(spring.flange_b, inertia3.flange_a)
                   connect(damper.flange_a, inertia2.flange_b)
                   connect(damper.flange_b, fixed.flange)
                   connect(sine.output, torque.tau)
                   connect(torque.support, fixed.flange)
                   connect(idealGear.support, fixed.flange)
                   connect(torque.flange, inertia1.flange_a)]

    @named model = ODESystem(connections, t,
        systems = [
            fixed,
            torque,
            inertia1,
            idealGear,
            inertia2,
            spring,
            inertia3,
            damper,
            sine
        ])
    @test_skip begin
        sys = structural_simplify(model) #key 7 not found
        prob = ODAEProblem(sys, Pair[], (0, 1.0))
        sol = solve(prob, Rodas4())
        @test SciMLBase.successful_retcode(sol)
    end
    # Plots.plot(sol; vars=[inertia2.w, inertia3.w])
end

@testset "Stick-Slip" begin
    @component function VelocityProfile(; name)
        @named sine = Blocks.Sine(amplitude = 10, frequency = 0.1)
        @named dz = Blocks.DeadZone(u_max = 2)
        @named lim = Blocks.Limiter(y_max = 6)
        @named output = Blocks.RealOutput()
        connections = [connect(sine.output, dz.input)
                       connect(dz.output, lim.input)
                       connect(lim.output, output)]
        ODESystem(connections, t, [], []; name = name, systems = [sine, dz, lim, output])
    end

    @named fixed = Fixed()
    @named spring = Spring(c = 6.5)
    @named damper = Damper(d = 0.01)
    @named inertia = Inertia(J = 0.0001)
    @named friction = RotationalFriction(f = 0.001, tau_c = 20, w_brk = 0.06035,
        tau_brk = 25)
    @named vel_profile = VelocityProfile()
    @named source = Speed()
    @named angle_sensor = AngleSensor()

    connections = [connect(vel_profile.output, source.w_ref)
                   connect(source.flange, friction.flange_a)
                   connect(friction.flange_b, inertia.flange_a)
                   connect(inertia.flange_b, spring.flange_a, damper.flange_a)
                   connect(spring.flange_b, damper.flange_b, fixed.flange)
                   connect(angle_sensor.flange, inertia.flange_a)]

    @named model = ODESystem(connections, t,
        systems = [
            fixed,
            inertia,
            spring,
            damper,
            vel_profile,
            source,
            friction,
            angle_sensor
        ])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, Pair[], (0, 10.0))

    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)
    @test sol[angle_sensor.phi.u] == sol[inertia.flange_a.phi]

    # p1 = Plots.plot(sol; vars=[inertia.flange_a.phi, source.phi], title="Angular Position", labels=["Inertia" "Source"], ylabel="Angle in rad")
    # p2 = Plots.plot(sol; vars=[friction.w_rel], title="Rel. Angular Velocity of Friction", label="", ylabel="Angular Velocity in rad/s")
    # Plots.plot(p1, p2, layout=(2, 1))
    # Plots.savefig("stick_slip.png")

    # Plots.scatter(sol[friction.w], sol[friction.tau], label="")
end

@testset "sensors" begin
    @named fixed = Fixed()
    @named inertia1 = Inertia(J = 2) # this one is fixed
    @named spring = Spring(c = 1e4)
    @named damper = Damper(d = 10)
    @named inertia2 = Inertia(J = 2, phi = pi / 2)
    @named speed_sensor = SpeedSensor()
    @named torque_sensor = TorqueSensor()
    @named rel_speed_sensor = RelSpeedSensor()

    connections = [connect(fixed.flange, inertia1.flange_b, rel_speed_sensor.flange_b)
                   connect(inertia1.flange_b, torque_sensor.flange_a)
                   connect(torque_sensor.flange_b, spring.flange_a, damper.flange_a,
                       speed_sensor.flange, rel_speed_sensor.flange_a)
                   connect(spring.flange_b, damper.flange_b, inertia2.flange_a)]

    @named model = ODESystem(connections,
        t,
        systems = [
            fixed, inertia1, inertia2, spring, damper, speed_sensor,
            rel_speed_sensor, torque_sensor
        ])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[], (0, 10.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[inertia1.w] .== 0)
    @test all(sol[inertia1.w] .== sol[speed_sensor.w.u])
    @test sol[inertia2.w][end]≈0 atol=1e-3 # all energy has dissipated
    @test all(sol[rel_speed_sensor.w_rel.u] .== sol[speed_sensor.w.u])
    @test all(sol[torque_sensor.tau.u] .== -sol[inertia1.flange_b.tau])

    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, Pair[], (0, 10.0))
    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[inertia1.w] .== 0)
    @test all(sol[inertia1.w] .== sol[speed_sensor.w.u])
    @test sol[inertia2.w][end]≈0 atol=1e-3 # all energy has dissipated
    @test all(sol[rel_speed_sensor.w_rel.u] .== sol[speed_sensor.w.u])
    @test all(sol[torque_sensor.tau.u] .== -sol[inertia1.flange_b.tau])

    # Plots.plot(sol; vars=[inertia1.w, inertia2.w])
end
