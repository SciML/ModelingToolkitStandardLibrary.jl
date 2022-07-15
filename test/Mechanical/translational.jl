using ModelingToolkitStandardLibrary.Mechanical.Translational, ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Blocks

@parameters t
D = Differential(t)

# see: https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Mechanics.Translational.Examples.Damper.html
@testset "mass damper" begin

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

    @named fixed = Fixed(; s0)
    @named mass = Mass(; m, s_start, v_start)
    @named damper = Damper(; d)

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
@testset "underdamped mass spring damper" begin

    mass_s_modelica = [0.100000000000000 3.93584712686420;
    0.860000000000000 3.91785659408340;
    1.04000000000000 2.79892156799871;
    1.46000000000000 1.47481226130478;
    1.96000000000000 2.80484384100597;
    2.46000000000000 3.92217116740851;
    2.92000000000000 3.22481919080783;
    3.46000000000000 2.44263894555251;
    3.94000000000000 2.88357609306897;
    4.46000000000000 3.33674208163243;
    4.94000000000000 3.07559491610384;
    5.46000000000000 2.79662733802557;
    5.92000000000000 2.94086200773602;
    6.42000000000000 3.12031467531611;
    7.44000000000000 2.92659078920977;
    8.48000000000000 3.04498080383203;
    9.40000000000000 2.97434929701123] # solver tolerance = 1e-12

    s0 = 4.5
    m = 1
    s_start = 3
    v_start = 10
    c = 10 # spring constant
    s_rel0 = 1.5 # unstretched spring length
    d = 1

    @named fixed = Fixed(; s0)
    @named mass = Mass(; m, s_start, v_start)
    @named damper = Damper(; d)
    @named spring = Spring(; c, s_rel0)
    
    connections = [
        connect(mass.flange_b, damper.flange_a, spring.flange_a)
        connect(damper.flange_b, spring.flange_b, fixed.flange)
    ]

    @named model = ODESystem(connections, t, systems=[fixed, mass, damper, spring])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 20.0), saveat=0.002)
    sol = solve(prob, DFBDF(), abstol=1e-12, reltol=1e-12)
    
    mass_s_julia = sol[mass.s]
    mass_v_julia = sol[mass.v]
    mass_s_julia_subset = mass_s_julia[in.(sol[t], [mass_s_modelica[:, 1]])]

    @test sol.retcode == :Success
    @test mass_v_julia[end] ≈ 0 atol = 1e-3 # all energy has dissipated
    @test mean(broadcast(abs, (mass_s_julia_subset - mass_s_modelica[:, 2]))) ≈ 0 atol = 1e-5 # jl and modelica s vs. t results are approx equal

end