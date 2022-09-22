
module TransPosExample
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
using Plots

@parameters t
D = Differential(t)

@named damping = Damper(d = 1)
@named spring = Spring(c = 1, s_rel0 = 0)
@named body = Mass(m = 1, v_start = 1, s_start = 3)
@named ground = Fixed(s0 = 3)

eqs = [connect(spring.flange_a, body.flange_a, damping.flange_a)
       connect(ground.flange, spring.flange_b, damping.flange_b)]

@named model = ODESystem(eqs, t; systems = [ground, body, spring, damping])

sys = structural_simplify(model)
vars = states(sys)
pars = parameters(sys)
du0 = D.(vars) .=> 0.0

prob = DAEProblem(sys, du0, [], (0, 10.0), [])
sol = solve(prob, DFBDF())

plot(sol, idxs = [body.s])

end
