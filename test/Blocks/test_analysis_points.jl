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
