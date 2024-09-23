using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

using ModelingToolkitStandardLibrary.Blocks: Parameter

@component function Step(;
        name, height = 1.0, offset = 0.0, start_time = 0.0, duration = Inf,
        smooth = 1e-5)
    @named output = B.RealOutput()
    duration_numeric = duration
    pars = @parameters offset=offset start_time=start_time height=height duration=duration step_val(t)::Bool=true
    equation = if smooth == false # use comparison in case smooth is a float
        offset +
        ifelse((step_val) & (t < start_time + duration), height, zero(height))
        #   ifelse((start_time <= t) & (t < start_time + duration), height, zero(height))
    else
        smooth === true && (smooth = 1e-5)
        if duration_numeric == Inf
            smooth_step(t, smooth, height, offset, start_time)
        else
            smooth_step(t, smooth, height, offset, start_time) -
            smooth_step(t, smooth, height, zero(start_time), start_time + duration)
        end
    end

    eqs = [
        output.u ~ equation
    ]

    compose(
        ODESystem(eqs,
            t,
            [],
            pars;
            name = name,
            continuous_events = [[t ~ start_time, t ~ start_time + duration] => (
                (i, u, p, c) -> i.ps[p.step_val] = false,
                [], [step_val], [step_val], nothing)]),
        [output])
end

#@testset "Fluid Domain and Tube" begin
function System(N; bulk_modulus, name)
    pars = @parameters begin
        bulk_modulus = bulk_modulus
    end

    systems = @named begin
        fluid = IC.HydraulicFluid(; bulk_modulus)
        stp = Step(;
            height = 2 * 101325, offset = 101325, start_time = 0.05, duration = Inf,
            smooth = false)
        src = IC.Pressure(;)
        vol = IC.FixedVolume(; vol = 10.0)
        res = IC.Tube(N; area = 0.01, length = 50.0)
    end

    eqs = [connect(stp.output, src.p)
           connect(fluid, src.port)
           connect(src.port, res.port_a)
           connect(res.port_b, vol.port)]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named sys5_1 = System(2; bulk_modulus = 1e9)

sys5_1 = structural_simplify(sys5_1)

initialization_eqs = [sys5_1.vol.port.p ~ 101325, sys5_1.res.p1.port_a.dm ~ 0]
initsys5_1 = ModelingToolkit.generate_initializesystem(sys5_1; initialization_eqs)

initsys5_1 = structural_simplify(initsys5_1)
initprob5_1 = NonlinearProblem(initsys5_1, [t => 0])
initsol5_1 = solve(initprob5_1)

prob5_1 = ODEProblem(sys5_1, [], (0, 50); initialization_eqs)
sol5_1 = solve(prob5_1) # ERROR: UndefRefError: access to undefined reference
