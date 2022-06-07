using ModelingToolkitStandardLibrary.Mechanical.Rotational, ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Blocks
# using Plots

@parameters t
D = Differential(t)

@testset "two inertias" begin
    @named fixed = Fixed()
    @named inertia1 = Inertia(J=2) # this one is fixed
    @named spring = Spring(c=1e4)
    @named damper = Damper(d=10)
    @named inertia2 = Inertia(J=2, phi_start=pi/2)

    connections = [
        connect(fixed.flange, inertia1.flange_b)
        connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
        connect(spring.flange_b, damper.flange_b, inertia2.flange_a)
    ]

    @named model = ODESystem(connections, t, systems=[fixed, inertia1, inertia2, spring, damper])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [], (0, 10.0))
    sol = solve(prob, DFBDF())

    # Plots.plot(sol; vars=[inertia1.w, inertia2.w])

    @test sol.retcode == :Success
    @test_skip all(sol[inertia1.w] .== 0)
    @test sol[inertia2.w][end] ≈ 0 atol=1e-3 # all energy has dissipated
end

@testset "two inertias with driving torque" begin
    amplitude = 10 # Amplitude of driving torque
    frequency = 5 # Frequency of driving torque
    J_motor = 0.1 # Motor inertia

    @named fixed = Fixed()
    @named torque = Torque(use_support=true)
    @named inertia1 = Inertia(J=2, phi_start=pi/2)
    @named spring = Spring(c=1e4)
    @named damper = Damper(d=10)
    @named inertia2 = Inertia(J=4)
    @named sine = Blocks.Sine(amplitude=amplitude, frequency=frequency)

    connections = [
        connect(sine.output, torque.tau)
        connect(torque.support, fixed.flange)
        connect(torque.flange, inertia1.flange_a)
        connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
        connect(spring.flange_b, damper.flange_b, inertia2.flange_a)
    ]

    @named model = ODESystem(connections, t, systems=[fixed, torque, inertia1, inertia2, spring, damper, sine])
    sys = structural_simplify(model)
    #prob = ODAEProblem(sys, Pair[], (0, 1.0))
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(inertia2.phi)) => 1.0], (0, 10.0))
    sol = solve(prob, DFBDF())

    # Plots.plot(sol; vars=[inertia1.w, -inertia2.w*2])

    @test_skip begin
        @test sol.retcode == :Success
        @test all(isapprox.(sol[inertia1.w], -sol[inertia2.w]*2, atol=1)) # exact opposite oscillation with smaller amplitude J2 = 2*J1
        @test_broken all(sol[torque.flange.tau] .== -sol[sine.output.u]) # torque source is equal to negative sine
    end
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
    @named torque = Torque(use_support=true)
    @named inertia1 = Inertia(J=J_motor)
    @named idealGear = IdealGear(ratio=ratio, use_support=true)
    @named inertia2 = Inertia(J=2)
    @named spring = Spring(c=1e4)
    @named inertia3 = Inertia(J=J_load)
    @named damper = Damper(d=damping)
    @named sine = Blocks.Sine(amplitude=amplitude, frequency=frequency)

    connections = [
        connect(inertia1.flange_b, idealGear.flange_a)
        connect(idealGear.flange_b, inertia2.flange_a)
        connect(inertia2.flange_b, spring.flange_a)
        connect(spring.flange_b, inertia3.flange_a)
        connect(damper.flange_a, inertia2.flange_b)
        connect(damper.flange_b, fixed.flange)
        connect(sine.output, torque.tau)
        connect(torque.support, fixed.flange)
        connect(idealGear.support, fixed.flange)
        connect(torque.flange, inertia1.flange_a)
    ]

    @named model = ODESystem(connections, t, systems=[fixed, torque, inertia1, idealGear, inertia2, spring, inertia3, damper, sine])
    sys = structural_simplify(model)
    @test_broken prob = ODAEProblem(sys, Pair[], (0, 1.0))
    # sol = solve(prob, Rodas4())
    # @test sol.retcode == :Success
    # Plots.plot(sol; vars=[inertia2.w, inertia3.w])
end
