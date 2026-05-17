using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Electrical.Machines
using SciCompDSL
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq: ReturnCode.Success
using OrdinaryDiffEqRosenbrock: Rodas4

# local mean (Statistics is not a test dependency of this package)
avg(x) = sum(x) / length(x)

# Machine under test (small surface PMSM)
R = 0.5
Ls = 1.0e-3
lambda = 0.01
pp = 4
Jr = 1.0e-4
Dr = 1.0e-4

@testset "PMSM: back-EMF constant" begin
    # Spin the shaft at a fixed speed; the PM back-EMF is then a balanced
    # 3-phase set of amplitude lambda*p*w0 (independent of stator load,
    # since e_x is an explicit function of speed and rotor angle).
    w0 = 100.0
    A = lambda * pp * w0                      # expected back-EMF amplitude

    @mtkmodel BackEMF begin
        @components begin
            machine = PMSM(
                R = R, Ls = Ls, lambda = lambda, p = pp,
                Jr = Jr, Dr = Dr
            )
            speed = Speed(exact = true)
            wref = Constant(k = w0)
            fixed = Fixed()
            ra = Resistor(R = 10.0)
            rb = Resistor(R = 10.0)
            rc = Resistor(R = 10.0)
            ground = Ground()
        end
        @equations begin
            connect(wref.output, speed.w_ref)
            connect(speed.flange, machine.flange)
            connect(machine.support, fixed.flange)
            connect(machine.pin_a, ra.p)
            connect(machine.pin_b, rb.p)
            connect(machine.pin_c, rc.p)
            connect(ra.n, rb.n, rc.n, ground.g)
        end
    end

    @mtkcompile sys = BackEMF()
    prob = ODEProblem(sys, unknowns(sys) .=> 0.0, (0.0, 0.1))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success

    # steady window (speed is imposed → reached immediately)
    mask = sol.t .> 0.03
    ea = sol[sys.machine.e_a][mask]
    eb = sol[sys.machine.e_b][mask]
    ec = sol[sys.machine.e_c][mask]

    # amplitude over a window spanning several electrical periods
    @test maximum(ea) ≈ A rtol = 2.0e-2
    @test minimum(ea) ≈ -A rtol = 2.0e-2
    # balanced 3-phase invariant (sampling-robust): Σ e² = 3/2 · A²
    @test avg(ea .^ 2 .+ eb .^ 2 .+ ec .^ 2) ≈ 1.5 * A^2 rtol = 2.0e-2
    # isolated wye: phase currents sum to zero
    sumi = sol[sys.machine.i_a] .+ sol[sys.machine.i_b] .+ sol[sys.machine.i_c]
    @test maximum(abs.(sumi)) < 1.0e-8
end

@testset "PMSM: locked-rotor q-axis torque" begin
    # Rotor locked at theta_e = 0; inject a balanced q-axis current set
    # (0, -Iq, +Iq). Closed form magnitude: |tau_e| = sqrt(3)·p·lambda·Iq.
    Iq = 5.0

    # Single line-to-line current source (pin_b -> pin_c) gives
    # i_b = +Iq, i_c = -Iq, i_a = 0 -- consistent with the isolated wye
    # (3 sources would make Σi = 0 redundant -> structurally singular).
    @mtkmodel LockedRotor begin
        @components begin
            machine = PMSM(
                R = R, Ls = Ls, lambda = lambda, p = pp,
                Jr = Jr, Dr = Dr
            )
            fixed_shaft = Fixed()
            fixed_house = Fixed()
            isrc = Current()
            iq = Constant(k = Iq)
            ground = Ground()
        end
        @equations begin
            connect(machine.flange, fixed_shaft.flange)
            connect(machine.support, fixed_house.flange)
            connect(iq.output, isrc.I)
            connect(isrc.p, machine.pin_b)
            connect(isrc.n, machine.pin_c)
            connect(machine.pin_a, ground.g)
        end
    end

    @mtkcompile sys = LockedRotor()
    prob = ODEProblem(sys, unknowns(sys) .=> 0.0, (0.0, 0.05))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success

    @test abs(sol[sys.machine.tau_e][end]) ≈ sqrt(3) * pp * lambda * Iq rtol = 1.0e-3
    # isolated wye: phase currents sum to zero
    @test abs(
        sol[sys.machine.i_a][end] + sol[sys.machine.i_b][end] +
            sol[sys.machine.i_c][end]
    ) < 1.0e-8
end
