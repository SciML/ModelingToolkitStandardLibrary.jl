using ModelingToolkitStandardLibrary.Magnetic, ModelingToolkit, OrdinaryDiffEq, Test
using SciCompDSL
using SciMLBase

import ModelingToolkitStandardLibrary.Electrical
import ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Magnetic
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq: ReturnCode.Success
using OrdinaryDiffEqRosenbrock: Rodas4
# using Plots

@testset "Inductor" begin
    mu_air = 1
    l_air = 0.0001
    mu_Fe = 1000
    l_Fe = 4 * 0.065
    a = b = 0.25

    @named source = Blocks.Sine(amplitude = 230 * sqrt(2), frequency = 50, phase = pi / 2)
    @named voltage = Electrical.Voltage()
    @named r = Electrical.Resistor(R = 7.5)
    @named ground = Electrical.Ground()
    @named coil = Magnetic.FluxTubes.ElectroMagneticConverter(N = 600, Phi = 0.0)
    @named ground_m = Magnetic.FluxTubes.Ground()
    @named r_mAirPar = Magnetic.FluxTubes.ConstantReluctance(R_m = a * b * l_air * mu_air)
    @named r_mFe = Magnetic.FluxTubes.ConstantReluctance(R_m = a * b * l_Fe * mu_Fe)
    @named r_mLeak = Magnetic.FluxTubes.ConstantReluctance(R_m = 1.2e6)
    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, r.p)
        connect(r.n, coil.p)
        connect(voltage.n, coil.n)
        connect(coil.port_p, r_mLeak.port_p)
        connect(r_mLeak.port_p, r_mAirPar.port_p)
        connect(r_mAirPar.port_n, r_mFe.port_p)
        connect(r_mFe.port_n, r_mLeak.port_n)
        connect(r_mFe.port_n, coil.port_n)
        connect(ground.g, voltage.n)
        connect(ground_m.port, r_mFe.port_n)
    ]
    @named model = System(
        connections, t,
        systems = [
            source,
            r,
            ground,
            coil,
            ground_m,
            r_mAirPar,
            r_mFe,
            r_mLeak,
            voltage,
        ]
    )
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0, 0.1))
    sol = solve(prob, Rodas4())

    # Plots.plot(sol; vars=[r.i])
    # Plots.plot(sol; vars=[r_mFe.V_m, r_mFe.Phi])

    @test SciMLBase.successful_retcode(sol)
    @test sol[r_mFe.Phi] == sol[r_mAirPar.Phi]
    coil_flux = sol[coil.port_p.Phi]
    leakage_flux = sol[r_mLeak.Phi]
    air_flux = sol[r_mAirPar.Phi]
    flux_balance = coil_flux + leakage_flux + air_flux
    flux_scale = maximum(maximum(abs, flux) for flux in (coil_flux, leakage_flux, air_flux))
    # Exact cancellation is lost when observables are reconstructed; four ULPs covers that
    # roundoff and the two additions.
    @test maximum(abs, flux_balance) <= 4 * eps(flux_scale)
end

@testset "EddyCurrent" begin
    # R is a dependent parameter bound to rho*l/A, so setting rho, l, A on a
    # problem must change the eddy-current dynamics.
    @named ec = Magnetic.FluxTubes.EddyCurrent(Phi = 0.0)
    @named src = Magnetic.FluxTubes.ConstantMagneticPotentialDifference(V_m = 1.0)
    @named g = Magnetic.FluxTubes.Ground()
    @named model = System(
        [
            connect(src.port_p, ec.port_p),
            connect(ec.port_n, src.port_n),
            connect(src.port_n, g.port),
        ], t; systems = [src, ec, g]
    )
    sys = mtkcompile(model)
    for (rho, l, A) in ((0.098e-6, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 2.0, 1 / 2))
        prob = ODEProblem(sys, [ec.rho => rho, ec.l => l, ec.A => A], (0.0, 1.0))
        @test prob.ps[ec.R] ≈ rho * l / A
        @test prob.f(prob.u0, prob.p, 0.0)[1] ≈ rho * l / A
    end
end
