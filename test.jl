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

N = 4
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
sol = solve(prob, Rodas5P());

import Plots
plt = Plots.plot();
for i in 2:N
    Plots.plot!(plt, sol, idxs=(links[i].x_2, links[i].y_2), lab=false)
end
lim = hypot(N, N)
Plots.plot!(plt, dpi=400, xlab="x", ylab="y", ylims=(-lim, lim), xlims=(-lim, lim))
Plots.savefig("$N-Pendulum.png")

using ModelingToolkitStandardLibrary.Blocks
@variables t
D = Differential(t)
@named fixed_y = TP.Fixed()
@named link = Link(; m = 1, l = 1, I = 1, g = -10)

@named mass = TP.Mass(m = 0.01, v0 = 0.0, s0 = 0)

eqs = [connect(link.TX1, mass.T)
       connect(link.TY1, fixed_y.T)]

@named model = ODESystem(eqs, t, [], []; systems = [link, fixed_y, mass])

sys = structural_simplify(model)
unset_diff = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
prob = ODEProblem(sys, [D(link.A) => 0; unset_diff .=> 0], (0, 10.0))
sol = solve(prob, Rodas5())
ts = range(0, 10, length=10 * 20)
us = sol(ts, idxs=[link.x_1, link.x_2, link.y_1, link.y_2]);
@gif for u in us
    plt = plot(u[1:2], u[3:4], lw = 2, lab=false, xlims=(-1.5, 0.5), ylims=(-1.5, 0.5), aspect_ratio=1, title="Pendulum linked to a mass")
end
