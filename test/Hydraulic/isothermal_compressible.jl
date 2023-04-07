using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

@parameters t
D = Differential(t)

# ------------------------------------------------------------
# Test fluid domain and Pipe
function System(N; bulk_modulus, name)
    pars = @parameters begin bulk_modulus = bulk_modulus end

    systems = @named begin
        fluid = IC.HydraulicFluid(; bulk_modulus)
        stp = B.Step(; height = 10e5, offset = 0, start_time = 0.005, duration = Inf,
                     smooth = true)
        src = IC.InputSource(; p_int = 0)
        vol = IC.FixedVolume(; p_int = 0, vol = 10.0)
    end

    if N == 1
        @named res = IC.PipeBase(; p_int = 0, area = 0.01, length = 500.0)
    else
        @named res = IC.Pipe(N; p_int = 0, area = 0.01, length = 500.0)
    end
    push!(systems, res)

    eqs = [connect(stp.output, src.input)
           connect(fluid, src.port)
           connect(src.port, res.port_a)
           connect(res.port_b, vol.port)]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named sys1_2 = System(1; bulk_modulus = 2e9)
@named sys1_1 = System(1; bulk_modulus = 1e9)
@named sys5_1 = System(5; bulk_modulus = 1e9)

NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 4 // 10)

syss = structural_simplify.([sys1_2, sys1_1, sys5_1])
probs = [ODEProblem(sys, [], (0, 0.2)) for sys in syss];
sols = [solve(prob, ImplicitEuler(nlsolve = NEWTON); initializealg = NoInit())
        for prob in probs];

s1_2 = complete(sys1_2)
s1_1 = complete(sys1_1)
s5_1 = complete(sys5_1)

# higher stiffness should compress more quickly and give a higher pressure
@test sols[1][s1_2.vol.port.p][end] > sols[2][s1_1.vol.port.p][end]

# N=5 pipe is compressible, will pressurize more slowly
@test sols[2][s1_1.vol.port.p][end] > sols[3][s5_1.vol.port.p][end]


# ------------------------------------------------------------
# Test the Valve 
function System(; name)
    pars = []

    systems = @named begin
        fluid = IC.HydraulicFluid()
        sink = IC.Source(; p = 10e5)
        vol = IC.FixedVolume(; vol = 0.1, p_int = 100e5)
        valve = IC.Valve(; p_a_int = 10e5, p_b_int = 100e5, area_int = 0, Cd = 1e6)
        ramp = B.Ramp(; height = 1, duration = 0.001, offset = 0, start_time = 0.001,
                      smooth = true)
    end

    eqs = [connect(fluid, sink.port)
           connect(sink.port, valve.port_a)
           connect(valve.port_b, vol.port)
           connect(valve.input, ramp.output)]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named valve_system = System()
s = complete(valve_system)
sys = structural_simplify(valve_system)
prob = ODEProblem(sys, [], (0, 0.01))
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt,
            initializealg = NoInit())

# the volume should discharge to 10bar
@test sol[s.vol.port.p][end] ≈ 10e5

# ------------------------------------------------------------
# Test DynmaicVolume and minimum_volume feature
function System(; name)
    pars = []

    # DynamicVolume values
    area = 0.01
    length = 0.1
    

    systems = @named begin
        fluid = IC.HydraulicFluid()
        src1 = IC.InputSource(; p_int = 10e5)
        src2 = IC.InputSource(; p_int = 10e5)

        vol1 = IC.DynamicVolume(; p_int = 10e5, area, direction = +1, dead_volume=2e-4, minimum_volume=2e-4)
        vol2 = IC.DynamicVolume(; p_int = 10e5, area, direction = -1, dead_volume=2e-4, minimum_volume=2e-4, x_int = length)
                                

        mass = T.Mass(; m = 100)

        sin1 = B.Sine(; frequency=0.5, amplitude = +1e5, offset=10e5)
        sin2 = B.Sine(; frequency=0.5, amplitude = -1e5, offset=10e5)

    end

    eqs = [connect(fluid, src1.port, src2.port)

           connect(src1.port, vol1.port)

           connect(src2.port, vol2.port)

           connect(vol1.flange, mass.flange, vol2.flange)
           
           connect(src1.input, sin1.output)
           connect(src2.input, sin2.output)
           ]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named min_vol_system = System()
s = complete(min_vol_system)
sys = structural_simplify(min_vol_system)
prob = ODEProblem(sys, [], (0, 5.))

sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt,
            initializealg = NoInit())

# begin
#     fig = Figure()

#     ax = Axis(fig[1,1], ylabel="position [m]", xlabel="time [s]")
#     lines!(ax, sol.t, sol[s.vol1.x]; label="vol1")
#     lines!(ax, sol.t, sol[s.vol2.x]; label="vol2")
    
#     ax = Axis(fig[2,1], ylabel="pressure [bar]", xlabel="time [s]")
#     lines!(ax, sol.t, sol[s.vol1.port.p]/1e5; label="vol1")
#     lines!(ax, sol.t, sol[s.vol2.port.p]/1e5; label="vol2")
#     Legend(fig[2,2], ax)
#     fig    
# end


# volume/mass should stop moving at opposite ends
@test round(sol[s.vol1.x][1]; digits = 2) == 0.0
@test round(sol[s.vol2.x][1]; digits = 2) == 0.1
@test round(sol[s.vol1.x][end]; digits = 2) == 0.1
@test round(sol[s.vol2.x][end]; digits = 2) == 0.0
