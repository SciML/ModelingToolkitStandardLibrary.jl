using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T




@parameters t
D = Differential(t)

function system(N; name, fluid)

    pars = @parameters begin
        fluid = fluid
    end

    systems = @named begin
        stp = B.Step(;height = 10e5, offset = 0, start_time = 0.005, duration = Inf, smooth = 0)
        src = IC.InputSource(;p_int=0)
        vol = IC.FixedVolume(;p_int=0, vol=10.0, fluid)
        
    end

    if N == 1
        @named res = IC.PipeBase(; p_int=0, area=0.01, length=500.0, fluid)
    else
        @named res = IC.Pipe(N; p_int=0, area=0.01, length=500.0, fluid)
    end
    push!(systems, res)


    eqs = [
        connect(stp.output, src.input)
        connect(src.port, res.port_a)
        connect(res.port_b, vol.port)
    ]

    ODESystem(eqs, t, [], pars; name, systems)
end


@named user_oil = IC.FluidType()
user_oil_type = typeof(IC.get_fluid(user_oil))
IC.density(::user_oil_type) = 876
IC.bulk_modulus(::user_oil_type) = 1.2e9
IC.viscosity(::user_oil_type) = 0.034


@named sys1 = system(1;fluid=IC.water_20C)
@named sys2 = system(2;fluid=IC.water_20C)
@named sys5 = system(5;fluid=IC.water_20C)
@named sys25 = system(25;fluid=IC.water_20C)


cc1 = structural_simplify(sys1; allow_parameter=false)
cc2 = structural_simplify(sys2; allow_parameter=false)
cc5 = structural_simplify(sys5; allow_parameter=false)
cc25 = structural_simplify(sys25; allow_parameter=false)

t_end = 0.2
prob0 = ODEProblem(cc1, [], (0,t_end), [complete(sys1).fluid => user_oil])
prob1 = ODEProblem(cc1, [], (0,t_end))
prob2 = ODEProblem(cc2, [], (0,t_end))
prob5 = ODEProblem(cc5, [], (0,t_end))
prob25 = ODEProblem(cc25, [], (0,t_end))


dt = 1e-4
NEWTON = NLNewton(check_div = false, always_new = true, max_iter=100, relax=4//10)
# sol0 = solve(prob0, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  
# sol1 = solve(prob1, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  
# sol2 = solve(prob2, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  
# sol5 = solve(prob5, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  
# sol25 = solve(prob25, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  




# using GLMakie
# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, sol0.t, sol0[sys1.res.port_a.p] .- sol0[sys1.res.port_b.p]; label="pipe base, oil", linestyle=:dash)
# lines!(ax, sol1.t, sol1[sys1.res.port_a.p] .- sol1[sys1.res.port_b.p]; label="pipe base")
# lines!(ax, sol2.t, sol2[sys2.res.port_a.p] .- sol2[sys2.res.port_b.p]; label="pipe N=2")
# lines!(ax, sol5.t, sol5[sys5.res.port_a.p] .- sol5[sys5.res.port_b.p]; label="pipe N=5")
# lines!(ax, sol25.t, sol25[sys25.res.port_a.p] .- sol25[sys25.res.port_b.p]; label="pipe N=25")
# axislegend(ax)
# fig




function actuator_system(; name, fluid_a, fluid_b = fluid_a)

    pars = @parameters begin
        fluid_a = fluid_a
        fluid_b = fluid_b
    end

    systems = @named begin
        stp = B.Step(;height = 10e5, offset = 0, start_time = 0.005, duration = Inf, smooth = 0)

        src_a = IC.InputSource(;p_int=0)
        res_a = IC.Pipe(5; p_int=0, area=0.01, length=500.0, fluid=fluid_a)
        act_a = IC.Actuator(+1; fluid=fluid_a, p_int=0, x_int=0, area=0.1, dead_volume=100)

        src_b = IC.InputSource(;p_int=0)
        res_b = IC.Pipe(5; p_int=0, area=0.01, length=500.0, fluid=fluid_b)
        act_b = IC.Actuator(-1; fluid=fluid_b, p_int=0, x_int=0, area=0.1, dead_volume=100)

        m = T.Mass(; m=1000)
        
    end




    eqs = [
        connect(stp.output, src_a.input, src_b.input)

        connect(src_a.port, res_a.port_a)
        connect(res_a.port_b, act_a.port)

        connect(src_b.port, res_b.port_a)
        connect(res_b.port_b, act_b.port)

        connect(act_a.flange, act_b.flange, m.flange)
    ]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named act_sys = actuator_system(; fluid_a = IC.water_20C)
cc = structural_simplify(act_sys; allow_parameter=false)

t_end = 0.5
prob1 = ODEProblem(cc, [], (0,t_end), [complete(act_sys).fluid_b => user_oil])
prob2 = ODEProblem(cc, [], (0,t_end))

sol1 = solve(prob1, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  
sol2 = solve(prob2, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit())  


using GLMakie
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, sol1.t, sol1[act_sys.act_a.x]; color=:red)
lines!(ax, sol1.t, sol1[act_sys.act_b.x]; color=:red, linestyle=:dash)
lines!(ax, sol2.t, sol2[act_sys.act_a.x]; color=:blue)
lines!(ax, sol2.t, sol2[act_sys.act_b.x]; color=:blue, linestyle=:dash)

ax = Axis(fig[2,1])
lines!(ax, sol1.t, (sol1[act_sys.act_a.port.p] .- sol1[act_sys.act_b.port.p])/1e5; color=:red)
lines!(ax, sol2.t, (sol2[act_sys.act_a.port.p] .- sol2[act_sys.act_b.port.p])/1e5; color=:blue)
fig