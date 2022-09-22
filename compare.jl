using ModelingToolkit
using OrdinaryDiffEq

@parameters t
D = Differential(t)

@connector function Port(; name, v0 = 0.0)
    pars = @parameters v0 = v0
    vars = @variables begin
        v(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, pars, name = name, defaults = Dict(v => v0, f => 0))
end

function Fixed(; name)
    @named port = Port()
    eqs = [port.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name), port)
end

function Body(; name, v0 = 0.0, m = 1.0, s0 = nothing, g = nothing)
    @named port = Port(v0 = v0) #TODO: accept parameters
    pars = @parameters m=m v0=v0
    vars = @variables begin
        v(t) = v0
        f(t) = 0
    end

    eqs = [port.v ~ v
           port.f ~ f]

    # gravity option
    if !isnothing(g)
        @parameters g = g
        push!(pars, g)
        push!(eqs, D(v) ~ f / m + g)
    else
        push!(eqs, D(v) ~ f / m)
    end

    # position option
    if !isnothing(s0)
        @parameters s0 = s0
        push!(pars, s0)

        @variables s(t) = s0
        push!(vars, s)

        push!(eqs, D(s) ~ v)
    end

    return compose(ODESystem(eqs, t, vars, pars; name = name), port)
end

function Spring(; name, k = 1e3, Δx0 = 0.0)

    #TODO: accept parameters
    @named port_a = Port()
    @named port_b = Port()

    pars = @parameters k = k
    vars = @variables begin
        Δx(t) = Δx0
        f(t) = k * Δx0
    end

    eqs = [D(Δx) ~ port_a.v - port_b.v
           f ~ k * Δx
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name), port_a, port_b)
end

function Damper(; name, d = 1e2, Δv0 = 0.0)

    #TODO: accept parameters
    @named port_a = Port()
    @named port_b = Port()

    pars = @parameters d = d
    vars = @variables v(t)=Δv0 f(t)=Δv0 * d

    eqs = [v ~ port_a.v - port_b.v
           f ~ v * d
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name), port_a, port_b)
end

@named damping = Damper(d = 1, Δv0 = 1)
@named spring = Spring(k = 1)
@named body = Body(m = 1, v0 = 1, s0 = 3)
@named ground = Fixed()

eqs = [connect(damping.port_a, body.port, spring.port_a)
       connect(ground.port, damping.port_b, spring.port_b)]

@named model = ODESystem(eqs, t; systems = [ground, body, damping, spring])

sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 10.0), [])
sol = solve(prob, ImplicitMidpoint(); dt = 0.01, initializealg = NoInit())

plot(sol, idxs = [body.s])

# defs = ModelingToolkit.defaults(sys)

# vars = states(model)

# for var in vars

#     x = sol[var][1]
#     x0 = defs[var]

#     while isa(x0, Sym)
#         x0 = defs[x0]
#     end

#     part1 = rpad("$var", 25)

#     println("$part1 $x \t $x0 \t def=$(defs[var])")

# end

var = spring.port_a.v
expr = ModelingToolkit.build_explicit_observed_function(sys, var; expression = true)

module TestPosition

using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition

@parameters t
D = Differential(t)

@named damping = Damper(d = 1)
@named spring = Spring(c = 1, s_rel0 = 3)
@named body = Mass(m = 1, v_start = 1, s_start = 3)
@named ground = Fixed()

eqs = [connect(damping.flange_a, body.flange_a, spring.flange_a)
       connect(ground.flange, damping.flange_b, spring.flange_b)]

@named model = ODESystem(eqs, t; systems = [ground, body, damping, spring])

sys = structural_simplify(model)
equations(sys)
prob = ODEProblem(sys, [], (0, 10.0), [])
sol = solve(prob, ImplicitMidpoint(); dt = 0.01, initializealg = NoInit())

plot(sol, idxs = [body.s])

end
