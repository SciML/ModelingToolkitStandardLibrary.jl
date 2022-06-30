import ModelingToolkitStandardLibrary.Blocks
include(pwd() * "/src/Mechanical/Translational/Translational.jl")
using .Translational, ModelingToolkit, OrdinaryDiffEq, Test
using CSV, DataFrames, Statisticss

@parameters t
D = Differential(t)

# see: https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Mechanics.Translational.Examples.Damper.html
@testset "fixed mass damper system" begin

    # Manually import select points from Modelica simulation (solver tolerance = 1e-12)

    mass_s_modelica = [0.002000000000000 3.019508803850819;
        0.010000000000000 3.088480745823828;
        0.022000000000000 3.169220984444292;
        0.038000000000000 3.245304115197263;
        0.066000000000000 3.323180201589973;
        0.092000000000000 3.359896450163389;
        0.108000000000000 3.373117745730875;
        0.124000000000000 3.381980258911502;
        0.142000000000000 3.388510086011018;
        0.164000000000000 3.393370882452241;
        0.182000000000000 3.395773080826793;
        0.218000000000000 3.398281442829992;
        0.382000000000000 3.399971513420815;
        0.788000000000000 3.400000000565588] # first column: time, second column: position

    # Initialize system state and define components and connections

    s0 = 4.5 # fixed offset position of flange housing
    m = 1 # mass of sliding mass
    s_start = 3 # initial value of absolute position of sliding mass
    v_start = 10 # initial value of absolute linear velocity of sliding mass
    d = 25 # damping constant

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
    sol = solve(prob, DFBDF(), abstol=1e-12, reltol=1e-12)
    
    mass_s_julia = sol[mass.s]
    mass_v_julia = sol[mass.v]
    mass_s_julia_subset = mass_s_julia[in.(sol[t], [mass_s_modelica[:, 1]])] # subset jl position array according to times defined in 'mass_s_modelica'

    # Define and run tests

    @test sol.retcode == :Success
    @test mass_v_julia[end] ≈ 0 atol = 1e-3 # all energy has dissipated
    @test mean(broadcast(abs, (mass_s_julia_subset - mass_s_modelica[:, 2]))) ≈ 0 atol = 1e-6 # jl and modelica s vs. t results are approx equal

end

# see: https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Mechanics.Translational.Examples.Damper.html
@testset "fixed mass spring damper system" begin

    mass_s_modelica = [0.008000000000000 3.072506925870173;
    0.016000000000000 3.131866370065987;
    0.028000000000000 3.201339780097094;
    0.048000000000000 3.279417481443546;
    0.092000000000000 3.359428741086432;
    0.112000000000000 3.374977582019366;
    0.130000000000000 3.383560768204081;
    0.144000000000000 3.387949651254135;
    0.162000000000000 3.391653220110489;
    0.194000000000000 3.395013356183882;
    0.346000000000000 3.395690582922517;
    1.800000000000000 3.373366033264920;
    9.400000000000000 3.275357787689842;
    26.600000000000000 3.138236634798095;
    54.600000000000000 3.045022807932026;
    89.600000000000000 3.011077566283765] # solver tolerance = 1e-12

    s0 = 4.5 # fixed offset position of flange housing
    m = 1 # mass of sliding mass
    s_start = 3 # initial value of absolute position of sliding mass
    v_start = 10 # initial value of absolute linear velocity of sliding mass
    d = 25 # damping constant
    c = 1 # spring constant
    s_rel0 = 1 # unstretched spring length

    @named fixed = Translational.Fixed(s0=s0)
    @named mass = Translational.Mass(m=m, s_start=s_start, v_start=v_start)
    @named damper = Translational.Damper(d=d)
    @named spring = Translational.Spring(c=c, s_rel0=s_rel0)

    connections = [
        connect(mass.flange_b, damper.flange_a)
        connect(mass.flange_b, spring.flange_a)
        connect(damper.flange_b, fixed.flange)
        connect(spring.flange_b, fixed.flange)
    ]

    @named model = ODESystem(connections, t, systems=[fixed, mass, damper, spring])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 180.0), saveat=0.002)
    sol = solve(prob, DFBDF(), abstol=1e-12, reltol=1e-12)
    
    mass_s_julia = sol[mass.s]
    mass_v_julia = sol[mass.v]
    mass_s_julia_subset = mass_s_julia[in.(sol[t], [mass_s_modelica[:, 1]])]

    @test sol.retcode == :Success
    @test mass_v_julia[end] ≈ 0 atol = 1e-3 # all energy has dissipated
    @test mean(broadcast(abs, (mass_s_julia_subset - mass_s_modelica[:, 2]))) ≈ 0 atol = 1e-6 # jl and modelica s vs. t results are approx equal

end