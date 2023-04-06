using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T
using IfElse

@parameters t
D = Differential(t)

NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 4 // 10)



function set_p_dm(;name, sgn=+1)
    
    pars = @parameters begin
        sgn = sgn
    end

    vars = @variables begin
        dm(t) = 0
    end
    
    systems = @named begin
        port = IC.HydraulicPort(; p_int=100e5)
    end

    eqs = [
        port.p ~ 100e5 + IfElse.ifelse(t<0.1, sgn*50e5*(t/0.1), sgn*50e5)
        port.dm ~ dm
    ]

    return ODESystem(eqs, t, vars, pars; name, systems)
end

function system_good(; name)
    pars = []
    vars = []

    systems = @named begin

        s1 = set_p_dm(; sgn = +1)
        s2 = set_p_dm(; sgn = -1)
        
        vol1 = IC.DynamicVolume(; p_int = 100e5, area = 0.01, direction = +1, dead_volume=0.01*0.1, minimum_volume=0.0001)
        vol2 = IC.DynamicVolume(; p_int = 100e5, area = 0.01, direction = -1, dead_volume=0.01*0.1, minimum_volume=0.0001)

        mass = T.Mass(; m=1000, s_0 = 0.0)
        
    end

    eqs = [
            connect(vol1.port, s1.port)
            connect(vol2.port, s2.port)

            connect(mass.flange, vol1.flange, vol2.flange)
           ]

    defaults = [
        s1.port.ρ => 1000
        s1.port.β => 1.2e9 
        s1.port.μ => 0.039
        s2.port.ρ => 1000
        s2.port.β => 1.2e9
        s2.port.μ => 0.039
        vol1.port.ρ => 1000
        vol1.port.β => 1.2e9
        vol1.port.μ => 0.039
        vol2.port.ρ => 1000
        vol2.port.β => 1.2e9
        vol2.port.μ => 0.039
    ]

    ODESystem(eqs, t, vars, pars; name, systems, defaults)
end

@named s = system_good()
s = complete(s)
sys = structural_simplify(s)
prob = ODEProblem(sys, [], (0., 0.5))
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); initializealg=NoInit(), adaptive=false, dt=1e-4)

# using CairoMakie
# lines(sol.t, sol[s.vol2.vol])


function system_bad(; name)
    pars = []
    vars = []

    systems = @named begin

        fluid = IC.HydraulicFluid()
        

        s1 = set_p_dm(; sgn = +1)
        s2 = set_p_dm(; sgn = -1)
        
        vol1 = IC.DynamicVolume(; p_int = 100e5, area = 0.01, direction = +1, dead_volume=0.01*0.1, minimum_volume=0.0001)
        vol2 = IC.DynamicVolume(; p_int = 100e5, area = 0.01, direction = -1, dead_volume=0.01*0.1, minimum_volume=0.0001)

        mass = T.Mass(; m=1000, s_0 = 0.0)
        
    end

    eqs = [
            
            connect(fluid, s1.port, s2.port)

            connect(vol1.port, s1.port)
            connect(vol2.port, s2.port)

            connect(mass.flange, vol1.flange, vol2.flange)
           ]



    ODESystem(eqs, t, vars, pars; name, systems)
end


@named s = system_bad()
s = complete(s)
sys = structural_simplify(s)
prob = ODEProblem(sys, [], (0., 0.5))
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); initializealg=NoInit(), adaptive=false, dt=1e-4)