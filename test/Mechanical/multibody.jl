using ModelingToolkit
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.Translational
using DifferentialEquations
# using Setfield

@parameters t

@named link = Link(;m=1, l=10, I=π*8.33, g=-9.807)
@named cart = Mass(;m=1, s_0=0)
# @named force = SineForce(;amp=3e3, freq=15)
# @named fixed = Fixed()
# @named m1 = Mass(;m=0.5)
# @named m2 = Mass(;m=0.5)

eqs = [
    connect(link.TX1, cart.flange) #, force.flange)
    # connect(link.TY1, fixed.flange)
    # connect(link.TX2, m1.flange)
    # connect(link.TY2, m2.flange)
]

@named model = ODESystem(eqs, t, [], []; systems=[link, cart])

sys = structural_simplify(model)
# sys = expand_connections(model)
# @set! sys.eqs = [(typeof(eq.lhs) != ModelingToolkit.Connection) & !ModelingToolkit._iszero(eq.lhs) & !ModelingToolkit.isdifferential(eq.lhs) ? 0 ~ eq.rhs - eq.lhs : eq for eq in sys.eqs]


prob = ODEProblem(sys, [], (0.0, 5), []; jac=true)
NEWTON = NLNewton(check_div = false, always_new = true, max_iter=100, relax=4//10)
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON), dt=0.0001, adaptive=false)
# sol = solve(prob)
plot(sol, idxs=[cart.s])