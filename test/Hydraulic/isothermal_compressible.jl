using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

@parameters t
D = Differential(t)

function system(N; bulk_modulus, name)
    pars = @parameters begin
        bulk_modulus = bulk_modulus
    end

    systems = @named begin
        fluid = IC.HydraulicFluid(;bulk_modulus)
        stp = B.Step(; height = 10e5, offset = 0, start_time = 0.005, duration = Inf, smooth = true)
        src = IC.InputSource(; p_int = 0)
        vol = IC.FixedVolume(; p_int = 0, vol = 10.0)
    end

    if N == 1
        @named res = IC.PipeBase(; p_int = 0, area = 0.01, length = 500.0)
    else
        @named res = IC.Pipe(N; p_int = 0, area = 0.01, length = 500.0)
    end
    push!(systems, res)

    eqs = [ connect(stp.output, src.input)
            connect(fluid, src.port)
            connect(src.port, res.port_a)
            connect(res.port_b, vol.port)]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named sys1_2 = system(1; bulk_modulus=2e9)
@named sys1_1 = system(1; bulk_modulus=1e9)
@named sys5_1 = system(5; bulk_modulus=1e9)

NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 4 // 10)

syss = structural_simplify.([sys1_2, sys1_1, sys5_1])
probs = [ODEProblem(sys, [], (0, 0.2)) for sys in syss];
sols = [solve(prob, ImplicitEuler(nlsolve=NEWTON); initializealg=NoInit()) for prob in probs];

s1_2 = complete(sys1_2)
s1_1 = complete(sys1_1)
s5_1 = complete(sys5_1)

# higher stiffness should compress more quickly and give a higher pressure
@test sols[1][s1_2.vol.port.p][end] > sols[2][s1_1.vol.port.p][end]

# N=5 pipe is compressible, will pressurize more slowly
@test sols[2][s1_1.vol.port.p][end] > sols[3][s5_1.vol.port.p][end]

# using CairoMakie
# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, sols[1].t, sols[1][s1_2.src.port.p]; label="sorce", color=:black)
# # lines!(ax, sols[2].t, sols[2][s1_1.src.port.p])
# scatterlines!(ax, sols[1].t, sols[1][s1_2.vol.port.p]; label="bulk=2e9")
# scatterlines!(ax, sols[2].t, sols[2][s1_1.vol.port.p]; label="bulk=1e9")
# scatterlines!(ax, sols[3].t, sols[3][s5_1.vol.port.p]; label="bulk=1e9, N=5")
# Legend(fig[1,2], ax)
# fig




function system(; bulk_modulus, name)
    pars = @parameters begin
        bulk_modulus = bulk_modulus
    end

    systems = @named begin
        fluid = IC.HydraulicFluid(;bulk_modulus)
        stp = B.Step(; height = 10e5, offset = 0, start_time = 0.05, duration = Inf, smooth = 0.001)
        src = IC.InputSource(; p_int = 0)
        vol1 = IC.DynamicVolume(; p_int = 0, area=0.01, direction=+1)
        vol2 = IC.DynamicVolume(; p_int = 0, area=0.01, direction=-1, x_int=1.0)
        mass = T.Mass(; m=1000, s_0=0)
        res = IC.PipeBase(; p_int = 0, area = 0.001, length = 50.0)
        cap = IC.Cap(;p_int=0)
    end

    eqs = [ connect(stp.output, src.input)
            connect(fluid, src.port, cap.port)
            connect(src.port, res.port_a)
            connect(res.port_b, vol1.port)
            connect(vol1.flange, mass.flange, vol2.flange)
            connect(vol2.port, cap.port)]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named sys1 = system(; bulk_modulus=1e9)
@named sys2 = system(; bulk_modulus=2e9)

syss = structural_simplify.([sys1, sys2])
probs = [ODEProblem(sys, [], (0, 0.2)) for sys in syss];
dt = 1e-4
sols = [solve(prob, ImplicitEuler(nlsolve=NEWTON); adaptive=false, dt, initializealg=NoInit()) for prob in probs];

s1 = complete(sys1)
s2 = complete(sys2)

# less stiff will compress more
@test maximum(sols[1][s1.mass.s]) > maximum(sols[2][s2.mass.s])


# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, sols[1].t, sols[1][s1.src.port.p]; label="sorce", color=:black)
# # lines!(ax, sols[2].t, sols[2][s1_1.src.port.p])
# scatterlines!(ax, sols[1].t, sols[1][s1.vol1.port.p]; label="bulk=1e9")
# scatterlines!(ax, sols[2].t, sols[2][s2.vol1.port.p]; label="bulk=2e9")
# Legend(fig[1,2], ax)
# fig


# fig = Figure()
# ax = Axis(fig[1,1])
# scatterlines!(ax, sols[1].t, sols[1][s1.mass.s]; label="bulk=1e9")
# scatterlines!(ax, sols[2].t, sols[2][s2.mass.s]; label="bulk=2e9")
# Legend(fig[1,2], ax)
# fig


# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, sols[1].t, sols[1][s1.vol1.dx]; label="bulk=2e9")
# lines!(ax, sols[1].t, sols[1][s1.vol2.dx]; label="bulk=2e9")
# lines!(ax, sols[1].t, sols[1][s1.mass.v]; label="bulk=2e9")
# Legend(fig[1,2], ax)
# fig
