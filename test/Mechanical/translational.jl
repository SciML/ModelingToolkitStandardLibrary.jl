import ModelingToolkitStandardLibrary.Blocks
include("/Users/kianmolani/Documents/Project Files/Codebases/ModelingToolkitStandardLibrary.jl/src/Mechanical/Translational/Translational.jl")
using ModelingToolkit, OrdinaryDiffEq, Test
using Plots

@parameters t
D = Differential(t)

function fixed_mass_damper_system()
    m = 1 # mass of sliding mass
    L = 1 # length of mass, from left flange to right flange (undefined)
    s_start = 3 # absolute position of sliding mass (ought to be of component center, if L were defined)
    v_start = 10 # absolute velocity of sliding mass
    d = 25 # damping constant of damper (expressed here in units of N*s/m)
    s0 = 4.5 # fixed offset position of housing

    @named fixed = Main.Translational.Fixed(s0=s0)
    @named mass = Main.Translational.Mass(m=m, s_start=s_start, v_start=v_start)
    @named damper = Main.Translational.Damper(d=d)

    connections = [
        connect(mass.flange_b, damper.flange_a)
        connect(damper.flange_b, fixed.flange)
    ]

    @named model = ODESystem(connections, t, systems=[fixed, mass, damper])
    sys = structural_simplify(model)
    prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 1.0))
    sol = solve(prob, DFBDF())

    display(Plots.plot(sol; vars=[mass.s]))
    display(Plots.plot(sol; vars=[mass.v]))
    display(Plots.plot(sol; vars=[mass.a]))
end

fixed_mass_damper_system()