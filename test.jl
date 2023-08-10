using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

dt = 4e-4
t_end = 10.0
time = 0:dt:t_end
x = @. time^2 + 1.0

@parameters t
D = Differential(t)

vars = @variables y(t)=1 dy(t)=0 ddy(t)=0
@named src = SampledData(; data=Float64[], dt)
@named int = Integrator()

eqs = [ y ~ src.output.u
        D(y) ~ dy
        D(dy) ~ ddy
        connect(src.output, int.input)]

@named sys = ODESystem(eqs, t, vars, []; systems = [int, src])
s = complete(sys)
sys = structural_simplify(sys)

prob = ODEProblem(sys, [], (0.0, t_end), [s.src.data => x]; tofloat=false)   
@time sol = solve(prob, ImplicitEuler());
# 0.001064 seconds (2.32 k allocations: 162.219 KiB)




# old SampledData using Parameter{Float64} type --------------------------------------------------
@named src = SampledData(Float64)
@named int = Integrator()

eqs = [ y ~ src.output.u
        D(y) ~ dy
        D(dy) ~ ddy
        connect(src.output, int.input)]

@named sys = ODESystem(eqs, t, vars, []; systems = [int, src])
s = complete(sys)
sys = structural_simplify(sys)

prob = ODEProblem(sys, [], (0.0, t_end), [s.src.buffer => Parameter(x,dt)]; tofloat=false)   
p = Parameter.(prob.p)
prob = remake(prob; p)
@time sol = solve(prob, ImplicitEuler());
# 0.000707 seconds (1.67 k allocations: 125.328 KiB)
