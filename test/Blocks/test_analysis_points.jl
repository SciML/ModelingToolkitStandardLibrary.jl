using Test, LinearAlgebra
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using ModelingToolkit: get_eqs, vars, @set!, get_iv

@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(-1)
t = ModelingToolkit.get_iv(P)

@test_logs (:warn,) (:warn,) connect(P.input, :bad_connection, C.output)

# Test with explicitly created AnalysisPoint
ap = AnalysisPoint(:plant_input)
eqs = [connect(P.output, C.input)
       connect(C.output, ap, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)

ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [P.x => 1], (0, 10))
sol = solve(prob, Rodas5())
@test norm(sol[1]) >= 1
@test norm(sol[end]) < 1e-6 # This fails without the feedback through C
# plot(sol)

matrices, _ = get_sensitivity(sys, ap)
@test matrices.A[] == -2
@test matrices.B[] * matrices.C[] == -1 # either one negative
@test matrices.D[] == 1

matrices, _ = get_comp_sensitivity(sys, ap)
@test matrices.A[] == -2
@test matrices.B[] * matrices.C[] == 1 # both positive or negative
@test matrices.D[] == 0

#=
# Equivalent code using ControlSystems. This can be used to verify the expected results tested for above.
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = 1                      # Negative feedback assumed in ControlSystems
S = sensitivity(P, C)      # or feedback(1, P*C)
T = comp_sensitivity(P, C) # or feedback(P*C)
=#

# Test with automatically created analysis point
eqs = [connect(P.output, C.input)
       connect(C.output, :plant_input, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)

matrices, _ = get_sensitivity(sys, :plant_input)
@test matrices.A[] == -2
@test matrices.B[] * matrices.C[] == -1 # either one negative
@test matrices.D[] == 1

matrices, _ = get_comp_sensitivity(sys, :plant_input)
@test matrices.A[] == -2
@test matrices.B[] * matrices.C[] == 1 # both positive
@test matrices.D[] == 0

## get_looptransfer

matrices, _ = Blocks.get_looptransfer(sys, :plant_input)
@test matrices.A[] == -1
@test matrices.B[] * matrices.C[] == -1 # either one negative
@test matrices.D[] == 0
#=
# Equivalent code using ControlSystems. This can be used to verify the expected results tested for above.
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = -1
L = P*C
=#

# Open loop
open_sys = Blocks.open_loop(sys, :plant_input)
@unpack u, y = open_sys

# Linearizing the open-loop system should yield the same system as get_looptransfer
matrices, _ = linearize(open_sys, [u], [y])
@test matrices.A[] == -1
@test matrices.B[] * matrices.C[] == -1 # either one negative
@test matrices.D[] == 0

# Test with more than one AnalysisPoint
eqs = [connect(P.output, :plant_output, C.input)
       connect(C.output, :plant_input, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)

matrices, _ = get_sensitivity(sys, :plant_input)
@test matrices.A[] == -2
@test matrices.B[] * matrices.C[] == -1 # either one negative
@test matrices.D[] == 1

## Test linearize between analysis points
matrices, _ = linearize(sys, :plant_input, :plant_output)
# Result should be the same as feedpack(P, 1), i.e., the closed-loop transfer function from plant input to plant output
@test matrices.A[] == -2
@test matrices.B[] * matrices.C[] == 1 # both positive
@test matrices.D[] == 0

## Test with subsystems

@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(1)
@named add = Blocks.Add(k2 = -1)
t = ModelingToolkit.get_iv(P)

eqs = [connect(P.output, :plant_output, add.input2)
       connect(add.output, C.input)
       connect(C.output, :plant_input, P.input)]

# eqs = [connect(P.output, add.input2)
#        connect(add.output, C.input)
#        connect(C.output, P.input)]

sys_inner = ODESystem(eqs, t, systems = [P, C, add], name = :inner)

@named r = Constant(k = 1)
@named F = FirstOrder(k = 1, T = 3)

eqs = [connect(r.output, F.input)
       connect(F.output, sys_inner.add.input1)]
sys_outer = ODESystem(eqs, t, systems = [F, sys_inner, r], name = :outer)

# test first that the structural_simplify works correctly
ssys = structural_simplify(sys_outer)
prob = ODEProblem(ssys, [P.x => 1], (0, 10))
# sol = solve(prob, Rodas5())
# plot(sol)

matrices, _ = get_sensitivity(sys_outer, :inner_plant_input)

using ControlSystemsBase # This is required to simplify the results to test against known solution
lsys = sminreal(ss(matrices...))
@test lsys.A[] == -2
@test lsys.B[] * lsys.C[] == -1 # either one negative
@test lsys.D[] == 1

matrices_So, _ = get_sensitivity(sys_outer, :inner_plant_output)
lsyso = sminreal(ss(matrices_So...))
@test lsys == lsyso || lsys == -1 * lsyso * (-1) # Output and input sensitivites are equal for SISO systems

## A more complicated test case
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks: t, Sine, PID, SecondOrder, Step, RealOutput
using ModelingToolkit: connect
# Parameters
m1 = 1
m2 = 1
k = 1000 # Spring stiffness
c = 10   # Damping coefficient
@named inertia1 = Inertia(; J = m1)
@named inertia2 = Inertia(; J = m2)
@named spring = Spring(; c = k)
@named damper = Damper(; d = c)
@named torque = Torque()

function SystemModel(u = nothing; name = :model)
    eqs = [connect(torque.flange, inertia1.flange_a)
           connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
           connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
    if u !== nothing
        push!(eqs, connect(torque.tau, u.output))
        return @named model = ODESystem(eqs, t;
                                        systems = [
                                            torque,
                                            inertia1,
                                            inertia2,
                                            spring,
                                            damper,
                                            u,
                                        ])
    end
    ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
end
function AngleSensor(; name)
    @named flange = Flange()
    @named phi = RealOutput()
    eqs = [phi.u ~ flange.phi
           flange.tau ~ 0]
    return ODESystem(eqs, t, [], []; name = name, systems = [flange, phi])
end

@named r = Step(start_time = 0)
model = SystemModel()
@named pid = PID(k = 100, Ti = 0.5, Td = 1)
@named filt = SecondOrder(d = 0.9, w = 10)
@named sensor = AngleSensor()
@named er = Add(k2 = -1)

connections = [connect(r.output, :r, filt.input)
               connect(filt.output, er.input1)
               connect(pid.ctr_output, :u, model.torque.tau)
               connect(model.inertia2.flange_b, sensor.flange)
               connect(sensor.phi, :y, er.input2)
               connect(er.output, :e, pid.err_input)]

closed_loop = ODESystem(connections, t, systems = [model, pid, filt, sensor, r, er],
                        name = :closed_loop)

prob = ODEProblem(structural_simplify(closed_loop), Pair[], (0.0, 4.0))
sol = solve(prob, Rodas4())
# plot(
#     plot(sol, vars = [filt.y, model.inertia1.phi, model.inertia2.phi]),
#     plot(sol, vars = [pid.ctr_output.u], title = "Control signal"),
#     legend = :bottomright,
# )

matrices, ssys = linearize(closed_loop, :r, :y)
lsys = ss(matrices...) |> sminreal
@test lsys.nx == 8

stepres = ControlSystemsBase.step(c2d(lsys, 0.001), 4)
@test stepres.y[:]≈sol(0:0.001:4, idxs = model.inertia2.phi) rtol=1e-4

# plot(stepres, plotx=true, ploty=true, size=(800, 1200), leftmargin=5Plots.mm)
# plot!(sol, vars = [model.inertia2.phi], sp=1, l=:dash)

matrices, ssys = get_sensitivity(closed_loop, :y)
So = ss(matrices...)

matrices, ssys = get_sensitivity(closed_loop, :u)
Si = ss(matrices...)

@test tf(So) ≈ tf(Si)
