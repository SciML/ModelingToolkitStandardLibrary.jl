module Example

using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkitStandardLibrary
using Plots

const TV = ModelingToolkitStandardLibrary.Mechanical.Translational
const TP = ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition

@parameters t
D = Differential(t)

# velocity --------------------------------------
@named damping = TV.Damper(d = 1)
@named spring1 = TV.Spring(k = 1, Δs0 = 0)
@named spring2 = TV.Spring(k = 2, Δs0 = 0)
@named body = TV.Mass(m = 1, v0 = 1, s0 = 3)
@named ground = TV.Fixed()

eqs = [connect(spring1.T1, body.T, damping.T1)
       connect(spring1.T2, damping.T2, spring2.T1)
       connect(spring2.T2, ground.T)]

@named model = ODESystem(eqs, t; systems = [ground, body, spring1, spring2, damping])

sys = structural_simplify(model)

#=
 Differential(t)(body₊v(t)) ~ (-damping₊d*body₊v(t) - spring₊k*spring₊Δs(t)) / body₊m
 Differential(t)(body₊s(t)) ~ body₊v(t)
 Differential(t)(spring₊Δs(t)) ~ body₊v(t)

 # m*a + d*v + k*Δs = 0
=#

prob = ODEProblem(sys, [], (0, 10.0), [])
sol = solve(prob, FBDF())

plot(sol, idxs = [body.s])

# position --------------------------------------
begin
    @named damping = TP.Damper(d = 1)
    @named spring1 = TP.Spring(k = 1, s1₀ = 3, s2₀ = 3)
    @named spring2 = TP.Spring(k = 2, s1₀ = 3, s2₀ = 3)
    @named body = TP.Mass(m = 1, v0 = 1, s0 = 3)
    @named ground = TP.Fixed(s0 = 3) # Note 1: Fixed position must be defined

    eqs = [connect(spring1.T1, body.T, damping.T1)
           connect(spring1.T2, damping.T2, spring2.T1)
           connect(ground.T, spring2.T2)]

    @named model = ODESystem(eqs, t; systems = [ground, body, spring1, spring2, damping])

    sys = structural_simplify(model)

    # Note 2: with ABS spring we get 3 equations, but based on Fixed position
    #=
    Differential(t)(body₊s(t)) ~ body₊v(t)
    Differential(t)(body₊v(t)) ~ (-ground₊T₊f(t)) / body₊m
    Differential(t)(body₊s(t)) ~ (ground₊T₊f(t) - spring₊k*(body₊s(t) - ground₊s0)) / damping₊d
    =#

    # Note 3: ERROR: The LHS operator must be unique. Please run `structural_simplify` or use the DAE form. Differential(t)(body₊s(t)) appears in LHS more than once.
    prob = ODEProblem(sys, [], (0, 10.0), [])

    # du0 = D.(states(sys)) .=> 0.0

    # Note 3: DAE must be used
    # prob = DAEProblem(sys, du0, [], (0, 10.0), [])
    sol = solve(prob, FBDF())

    plot(sol, idxs = [body.s])
end

# position --------------------------------------
begin
    @named damping = TP.Damper(d = 1)
    @named spring = TP.Spring(TP.REL; k = 1, Δs0 = 0)
    @named body = TP.Mass(m = 1, v0 = 1, s0 = 3)
    @named ground = TP.Fixed() # Note 4: using relative spring removes fixed position from equations

    eqs = [connect(spring.T1, body.T, damping.T1)
           connect(ground.T, spring.T2, damping.T2)]

    @named model = ODESystem(eqs, t; systems = [ground, body, spring, damping])

    sys = structural_simplify(model)

    # Note 5: relative spring adds additional equation
    #=

     Differential(t)(body₊s(t)) ~ body₊v(t)
     Differential(t)(body₊v(t)) ~ (-ground₊T₊f(t)) / body₊m
     Differential(t)(spring₊Δs(t)) ~ (ground₊T₊f(t) - spring₊k*spring₊Δs(t)) / damping₊d
     Differential(t)(body₊s(t)) ~ (ground₊T₊f(t) - spring₊k*spring₊Δs(t)) / damping₊d

    =#

    # du0 = D.(states(sys)) .=> 0.0
    prob = ODEProblem(sys, [], (0, 10.0), [])
    sol = solve(prob, FBDF())

    plot(sol, idxs = [body.s])
end
end
