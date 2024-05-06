using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using OrdinaryDiffEq: ReturnCode.Success
using Test

#=
Testing strategy:
The general strategy is to test systems using simple inputs where the solution
is known on closed form. For algebraic systems (without differential variables),
an integrator with a constant input is often used together with the system under test.
=#

@testset "Integrator" begin
    clock = Clock(t, 0.1)
    k = ShiftIndex(clock)
    @named c = Constant(; k = 1)
    @named sampler = Sampler(; clock)
    @named intc = Integrator()

    # Backward
    @named int = DiscreteIntegrator(x = 1)
    @named iosys = ODESystem(
        [connect(c.output, sampler.input)
         connect(sampler.output, int.input)
         connect(c.output, intc.input)],
        t,
        systems = [sampler, int, intc, c])
    model = complete(iosys)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[int.x(k - 1) => 0
                                int.u(k - 1) => 0], (0.0, 1.0))
    sol = solve(prob, Rodas4(), kwargshandle = KeywordArgSilent)
    @test sol.retcode == Success
    @test sol.prob.kwargs[:disc_saved_values][1].t ≈ 0:sampletime(clock):1
    @test reduce(vcat, sol.prob.kwargs[:disc_saved_values][1].saveval) ≈
          range(1.1, step = sampletime(clock), length = 11)

    # Forward
    @named int = DiscreteIntegrator(x = 1, method = :forward)
    @named iosys = ODESystem(
        [connect(c.output, sampler.input)
         connect(sampler.output, int.input)
         connect(c.output, intc.input)],
        t,
        systems = [sampler, int, intc, c])
    model = complete(iosys)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[int.x(k - 1) => 0
                                int.u(k - 1) => 0], (0.0, 1.0))
    sol = solve(prob, Rodas4(), kwargshandle = KeywordArgSilent)
    @test sol.retcode == Success
    @test sol.prob.kwargs[:disc_saved_values][1].t ≈ 0:sampletime(clock):1
    @test reduce(hcat, sol.prob.kwargs[:disc_saved_values][1].saveval)'[:, 2] ≈
          range(1.0, step = sampletime(clock), length = 11)

    # Tustin
    @named int = DiscreteIntegrator(x = 1, method = :tustin)
    @named iosys = ODESystem(
        [connect(c.output, sampler.input)
         connect(sampler.output, int.input)
         connect(c.output, intc.input)],
        t,
        systems = [sampler, int, intc, c])
    model = complete(iosys)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[int.x(k - 1) => 0
                                int.u(k - 1) => 0], (0.0, 1.0))
    sol = solve(prob, Rodas4(), kwargshandle = KeywordArgSilent)
    @test sol.retcode == Success
    @test sol.prob.kwargs[:disc_saved_values][1].t ≈ 0:sampletime(clock):1
    @test reduce(hcat, sol.prob.kwargs[:disc_saved_values][1].saveval)'[:, 2] ≈
          range(1.05, step = sampletime(clock), length = 11)
end

# for v in 𝑑vertices(graph)
#     vd = var_domain[v]
#     eq_inds = 𝑑neighbors(graph, v)
#     isempty(eq_inds) && continue
#     # eqs = collect(ts.sys.eqs)
#     for (i, eqi) in enumerate(eq_inds)
#         eq_domain[eqi] = vd

#         #     eq = eqs[eqi]
#         #     sampletime_vars = OrderedSet()
#         #     vars!(sampletime_vars, eq; op = InferredSampleTime2)
#         #     Ts = sampletime(eq_domain[eqi])
#         #     @show eqi, eq_domain[eq_inds[i]], Ts
#         #     @show substitute(eq, Dict(InferredSampleTime2() => Ts))
#         #     # @set! ts.sys.eqs = eqs[ieqs]
#         #     eqs[eqi] = substitute(eq, Dict(InferredSampleTime2() => Ts))
#     end
#     # @set! ts.sys.eqs = eqs
# end

# ==============================================================================
## DiscretePIDParallel
# ==============================================================================

dt = 0.05
k = ShiftIndex()

@mtkmodel ClosedLoop begin
    @components begin
        plant = FirstOrder(k = 1, T = 1)
        sampler = Blocks.Sampler(; dt)
        zoh = Blocks.ZeroOrderHold()
        controller = Blocks.DiscretePIDParallel(
            kp = 2, ki = 2, Imethod = :forward, with_D = false)
        ref = Constant(k = 0.5)
    end
    @equations begin
        connect(ref.output, controller.reference)
        connect(controller.ctr_output, zoh.input)
        connect(zoh.output, plant.input)
        connect(plant.output, sampler.input)
        connect(sampler.output, controller.measurement)
    end
end

#
@named model = ClosedLoop()
model = complete(model)
# ci, varmap = infer_clocks(expand_connections(model))
model = structural_simplify(model)

Tf = 5
timevec = 0:(dt):Tf

import ControlSystemsBase as CS
import ControlSystemsBase: c2d, tf, feedback, lsim
P = CS.c2d(CS.ss([-1], [1], [1], 0), dt)
C = CS.c2d(CS.ss([0], [1], [2], [2]), dt, :fwdeuler)

# Test the output of the continuous partition
G = feedback(P * C)
res = lsim(G, (x, t) -> [0.5], timevec)
y = res.y[:]

prob = ODEProblem(model,
    [model.plant.x => 0.0; model.controller.kp => 2.0; model.controller.ki => 2.0;
     model.controller.eI => 0.0],
    (0.0, Tf))

sol = solve(prob,
    Tsit5(),
    kwargshandle = KeywordArgSilent,
    abstol = 1e-8,
    reltol = 1e-8)

# plot(timevec, [y sol(timevec, idxs = model.plant.output.u)[:]], m = :o, lab = ["CS" "MTK"])
# display(current())

@test sol(timevec, idxs = model.plant.output.u)[:]≈y rtol=1e-6

@test_skip begin
    # Test the output of the discrete partition
    G = feedback(C, P)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]
    @test_broken sol(timevec .+ 1e-10, idxs = model.controller.output.u)≈y rtol=1e-8 # Broken due to discrete observed
    # plot([y sol(timevec .+ 1e-12, idxs=model.controller.output.u)], lab=["CS" "MTK"])
end

# ==============================================================================
## DiscretePIDStandard
# ==============================================================================

dt = 0.05
k = ShiftIndex()

@mtkmodel ClosedLoop begin
    @components begin
        plant = FirstOrder(k = 1, T = 1)
        sampler = Blocks.Sampler(; dt)
        zoh = Blocks.ZeroOrderHold()
        controller = Blocks.DiscretePIDStandard(
            K = 2, Ti = 1, Imethod = :forward, with_D = false)
        ref = Constant(k = 0.5)
    end
    @equations begin
        connect(ref.output, controller.reference)
        connect(controller.ctr_output, zoh.input)
        connect(zoh.output, plant.input)
        connect(plant.output, sampler.input)
        connect(sampler.output, controller.measurement)
    end
end

#
@named model = ClosedLoop()
model = complete(model)
# ci, varmap = infer_clocks(expand_connections(model))
model = structural_simplify(model)

Tf = 5
timevec = 0:(dt):Tf

import ControlSystemsBase as CS
import ControlSystemsBase: c2d, tf, feedback, lsim
P = CS.c2d(CS.ss([-1], [1], [1], 0), dt)
C = CS.c2d(CS.ss([0], [1], [2], [2]), dt, :fwdeuler)

# Test the output of the continuous partition
G = feedback(P * C)
res = lsim(G, (x, t) -> [0.5], timevec)
y = res.y[:]

prob = ODEProblem(model,
    [model.plant.x => 0.0; model.controller.eI => 0.0],
    (0.0, Tf))

sol = solve(prob,
    Tsit5(),
    kwargshandle = KeywordArgSilent,
    abstol = 1e-8,
    reltol = 1e-8)

# plot(timevec, [y sol(timevec, idxs = model.plant.output.u)[:]], m = :o, lab = ["CS" "MTK"])
# display(current())

@test sol(timevec, idxs = model.plant.output.u)[:]≈y rtol=1e-6

@test_skip begin
    # Test the output of the discrete partition
    G = feedback(C, P)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]
    @test_broken sol(timevec .+ 1e-10, idxs = model.controller.output.u)≈y rtol=1e-8 # Broken due to discrete observed
    # plot([y sol(timevec .+ 1e-12, idxs=model.controller.output.u)], lab=["CS" "MTK"])
end

# ==============================================================================
## Delay
# ==============================================================================

@mtkmodel DelayModel begin
    @components begin
        fake_plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
        input = Step(start_time = 2, smooth = false)
        sampler = Sampler(; dt = 1)
        delay = Delay(n = 3)
        zoh = ZeroOrderHold()
    end
    @equations begin
        connect(input.output, sampler.input)
        connect(sampler.output, delay.input)
        connect(delay.output, zoh.input)
        connect(zoh.output, fake_plant.input)
    end
end

@mtkbuild m = DelayModel()
prob = ODEProblem(
    m, [m.delay.u(k - 3) => 0, m.delay.u(k - 2) => 0, m.delay.u(k - 1) => 0], (0.0, 10.0))
sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)

@test reduce(vcat, sol((0:10) .+ 1e-2))[:]≈[zeros(5); ones(6)] atol=1e-2

# ==============================================================================
## Difference
# ==============================================================================
using ModelingToolkitStandardLibrary.Blocks
k = ShiftIndex(Clock(t, 1))

@mtkmodel DiffModel begin
    @components begin
        input = Step(start_time = 2, smooth = false)
        diff = Blocks.Difference(z = k)
        zoh = Blocks.ZeroOrderHold()
        plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
    end
    @equations begin
        connect(input.output, diff.input)
        connect(diff.output, zoh.input)
        connect(zoh.output, plant.input)
    end
end

@mtkbuild m = DiffModel()
prob = ODEProblem(m, Dict(m.diff.u(k - 1) => 0), (0.0, 10.0))
sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)
@test reduce(vcat, sol((0:10) .+ 1e-2))[:]≈[zeros(2); 1; zeros(8)] atol=1e-2
