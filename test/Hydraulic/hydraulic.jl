using ModelingToolkitStandardLibrary.Hydraulic
using ModelingToolkit, OrdinaryDiffEq, Test
# using Plots

@parameters t
D = Differential(t)

@testset "restriction and volume (shows dynamic compressibility)" begin
    @named reservoir_left = Reservoir(; p = 2e5)
    @named valve = LocalRestriction(; A = 1e-4, Cd = 0.64)
    @named volume = ConstantVolume(; V = 0.001)
    connections = [
        connect(reservoir_left.a, valve.a)
        connect(valve.b, volume.a)
    ]
    @named model = ODESystem(connections, t; systems = [reservoir_left, valve, volume])

    # Water
    @named fluid_props = FluidProperties(;
        p_atm = 101325.0,
        nu_atm = 1.0034e-6,
        beta_atm = 2.1791e9,
        rho_atm = 998.21,
    )
    model = extend(model, fluid_props)

    sys = structural_simplify(model)
    prob = ODEProblem(sys, [volume.p => 1e5], (0.0, 0.5e-3))
    sol = solve(prob, Tsit5(), reltol = 1e-6)

    @test sol.retcode == :Success
    @test sol[volume.p][end] ≈ 2e5 atol = 1
    @test sol[volume.a.m_flow][end] ≈ 0 atol = 1e-2

    # p1=plot(sol; vars=[volume.p / 1e6], ylabel="p in MPa", label="Pressure of Volume")
    # p2=plot(sol; vars=[volume.a.m_flow], ylabel="m_flow in kg/s", label="Mass Flow Rate into the volume")
    # plot(p1, p2, layout=(2,1))
end