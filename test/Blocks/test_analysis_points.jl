using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using ModelingToolkit: get_eqs, vars, @set!, get_iv

@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(-1)
t = ModelingToolkit.get_iv(P)

# Test with explicitly created AnalysisPoint
ap = AnalysisPoint(:plant_input)
eqs = [connect(P.output, C.input)
       connect(C.output, ap, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)

ssys = structural_simplify(expand_analysis_points(sys))
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
@test matrices.A[] == -1
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

# test first that the structural_simplify ∘ expand_analysis_points works correctly
ssys = structural_simplify(expand_analysis_points(sys_outer))
# ssys = structural_simplify((sys_outer))
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
