using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B

@parameters t
D = Differential(t)

function system(; name, fluid)

    pars = @parameters begin
        fluid = fluid
    end

    systems = @named begin
        stp = B.Step(;height = 100e5, offset = 0, start_time = 0.01, duration = Inf, smooth = 1e-5)
        src = IC.InputSource(;p=0)
        vol = IC.FixedVolume(;p_int=0, vol=0.1, fluid)
        res = IC.LaminarResistance(; p_int=0, area=0.01, length=1e2, fluid)
    end

    eqs = [
        connect(stp.output, src.input)
        connect(src.port, res.port_a)
        connect(res.port_b, vol.port)
    ]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named sys = system(;fluid=IC.Fluids.water)

cc = structural_simplify(sys; allow_parameter=false)
prob = ODEProblem(cc, [], (0,0.1))
dt = 1e-6
sol = solve(prob, ImplicitEuler(); adaptive=false, dt, initializealg=NoInit())  

# sys = complete(sys)
# sol[sys.vol.port.p] 

# using GLMakie

# lines(sol.t, sol[sys.vol.port.p]/1e5)

# function ModelingToolkit.promote_to_concrete(vs; tofloat = true, use_union = false)
#     if isempty(vs)
#         return vs
#     end
#     T = eltype(vs)
#     if Base.isconcretetype(T) && (!tofloat || T === float(T)) # nothing to do
#         vs
#     else
#         sym_vs = filter(x -> SymbolicUtils.issym(x) || SymbolicUtils.istree(x), vs)
#         isempty(sym_vs) || throw_missingvars_in_sys(sym_vs)
#         C = typeof(first(vs))
#         println("55: C=$C")
#         I = Int8
#         has_int = false
#         has_array = false
#         array_T = nothing
#         for v in vs
#             if v isa AbstractArray
#                 has_array = true
#                 array_T = typeof(v)
#             end
#             E = eltype(v)
#             print("67: C (before)=$C")
#             C = promote_type(C, E)
#             println(", C (after)=$C, E=$E")
#             if E <: Integer
#                 has_int = true
#                 I = promote_type(I, E)
#             end
#         end
#         if tofloat && !has_array
#             C = float(C)
#         elseif has_array || (use_union && has_int && C !== I)
#             if has_array
#                 C = Union{C, array_T}
#             end
#             if has_int
#                 C = Union{C, I}
#             end
#             return copyto!(similar(vs, C), vs)
#         end
#         convert.(C, vs)
#     end
# end

# ModelingToolkit.promote_to_concrete(prob.p)