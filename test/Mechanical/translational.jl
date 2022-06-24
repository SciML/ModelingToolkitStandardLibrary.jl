import ModelingToolkitStandardLibrary.Blocks
include(pwd() * "/src/Mechanical/Translational/Translational.jl")
using .Translational, ModelingToolkit, OrdinaryDiffEq, Test
using CSV, DataFrames, Statistics

@parameters t
D = Differential(t)

# see: https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Mechanics.Translational.Examples.Damper.html
@testset "fixed mass damper system" begin

    # Import Modelica simulation data

    df = CSV.read("./test/Mechanical/modelica/FixedMassDamper_Linear.csv", DataFrame)
    df = df[Not([502]), :]

    time = df.time
    mass_s_modelica = df."mass1.s"
    mass_v_modelica = df."mass1.v"
    mass_a_modelica = df."mass1.a"

    # Initialize system state and define components and connections

    s0 = 4.5 # fixed offset position of flange housing
    m = 1 # mass of sliding mass
    s_start = 3 # initial value of absolute position of sliding mass
    v_start = 10 # initial value of absolute linear velocity of sliding mass
    d = 25 # damping constant of damper

    @named fixed = Translational.Fixed(s0=s0)
    @named mass = Translational.Mass(m=m, s_start=s_start, v_start=v_start)
    @named damper = Translational.Damper(d=d)

    connections = [
        connect(mass.flange_b, damper.flange_a)
        connect(damper.flange_b, fixed.flange)
    ]

    # Solve ODE system

    @named model = ODESystem(connections, t, systems=[fixed, mass, damper])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 1.0), saveat=0.002)
    sol = solve(prob, DFBDF())

    mass_s_julia = sol[mass.s]
    mass_v_julia = sol[mass.v]
    mass_a_julia = sol[mass.a]

    # Define and run tests

    @test sol.retcode == :Success
    @test mass_v_julia[end] ≈ 0 atol = 1e-3 # all energy has dissipated
    @test mean(broadcast(abs, (mass_s_julia - mass_s_modelica))) ≈ 0 atol = 1e-5 # jl and modelica s vs. t results are approx equal
    @test mean(broadcast(abs, (mass_v_julia - mass_v_modelica))) ≈ 0 atol = 1e-4 # jl and modelica v vs. t results are approx equal
    @test mean(broadcast(abs, (mass_a_julia - mass_a_modelica))) ≈ 0 atol = 1e-2 # jl and modelica a vs. t results are approx equal
end

