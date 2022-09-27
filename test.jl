using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.Translational: Fixed
const TP = ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
using ModelingToolkit
using OrdinaryDiffEq

@variables t
D = Differential(t)
@named fixed1 = TP.Fixed()
@named fixed2 = TP.Fixed()
@named link = Link(; m = 1, l = 1, I = 1, g = -10)

eqs = [connect(link.TX1, fixed1.T)
       connect(link.TY1, fixed2.T)]

@named model = ODESystem(eqs, t, [], []; systems = [link, fixed1, fixed2])

sys = structural_simplify(model)
prob = ODEProblem(sys, [D(link.A) => 0], (0, 10.0))
sol = solve(prob, Rodas5())

N = 20
@named links[1:N] = Link(; m = 1, l = 1, I = 1, g = -10);

eqs = [
       connect(links[1].TX1, fixed1.T)
       connect(links[1].TY1, fixed2.T)
      ]
for i in 2:N
    push!(eqs, connect(links[i].TX1, links[i-1].TX2))
    push!(eqs, connect(links[i].TY1, links[i-1].TY2))
end

@named model = ODESystem(eqs, t, [], []; systems = [links; fixed1; fixed2]);

sys = structural_simplify(model)
unset_diff = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
u0 = [s => pi/i for (i, s) in enumerate(states(sys))]
prob = ODEProblem(sys, [unset_diff .=> 0.0; u0], (0, 10.0));
#sol = solve(prob, Rodas5P());
sol = solve(prob, ImplicitEuler(nlsolve = NLNewton(check_div = false, always_new = true)), dt=0.001, reltol=1e-5, adaptive=false);

import Plots
plt = Plots.plot();
for i in 2:N
    Plots.plot!(plt, sol, idxs=(links[i].x_2, links[i].y_2), lab=false)
end
lim = hypot(N, N)
Plots.plot!(plt, dpi=400, xlab="x", ylab="y", ylims=(-lim, lim), xlims=(-lim, lim))
Plots.savefig("$N-Pendulum.png")
