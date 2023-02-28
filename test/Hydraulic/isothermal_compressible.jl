using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC

@parameters t
D = Differential(t)

function system(; name, fluid)

    pars = @parameters begin
        p = 100e5
    end

    systems = @named begin
        src = IC.Source(;p)
        vol = IC.FixedVolume(;p_int=p, vol=0.1, fluid)
        res = IC.LaminarResistance(; p_int=p, area=0.01, length=0.1, fluid)
    end

    eqs = [
        connect(src.port, res.port_a)
        connect(res.port_b, vol.port)
    ]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named sys = system(;fluid=Int(IC.water))
ss = structural_simplify(sys)

# prob = ODEProblem(ss, [], (0,0.01))
# sol = solve(prob, ImplicitEuler())