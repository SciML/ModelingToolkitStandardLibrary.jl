using ModelingToolkit: get_eqs, vars, @set!, get_iv

Base.@kwdef mutable struct AnalysisPoint
    in = nothing
    out = nothing
    name::Symbol
end

"""
    AnalysisPoint(in, out, name::Symbol)
    AnalysisPoint(in, out; name::Symbol)
    AnalysisPoint(name::Symbol)

Create an AnalysisPoint for linear analysis. Analysis points can also be created automatically by calling 
```
connect(in, :ap_name, out)
```

# Arguments:
- `in`: A connector of type [`RealOutput`](@ref).
- `out`: A connector of type [`RealInput`](@ref).
- `name`: The name of the analysis point.

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref), [`open_loop`](@ref)

# Example
```julia
using ModelingToolkitStandardLibrary.Blocks
@named P = FirstOrder(k=1, T=1)
@named C = Gain(-1)
t = ModelingToolkit.get_iv(P)
eqs = [
    connect(P.output, C.input)
    connect(C.output, :plant_input, P.input)   # Connect with an automatically created analysis point
]
sys = ODESystem(eqs, t, systems=[P,C], name=:feedback_system)

matrices_S, _ = get_sensitivity(sys, :plant_input) # Compute the matrices of a state-space representation of the (input) sensitivity funciton.
matrices_T, _ = get_comp_sensitivity(sys, :plant_input)
```
Continued linear analysis and design can be performed using ControlSystemsBase.jl.
Create `ControlSystemsBase.StateSpace` objects using
```julia
using ControlSystemsBase, Plots
S = ss(matrices_S...)
T = ss(matrices_T...)
bodeplot([S, T], lab=["S" "T"])
```
The sensitivity functions obtained this way should be equivalent to the ones obtained with the code below
```
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = 1                      # Negative feedback assumed in ControlSystems
S = sensitivity(P, C)      # or feedback(1, P*C)
T = comp_sensitivity(P, C) # or feedback(P*C)
```
"""
AnalysisPoint(in, out; name) = AnalysisPoint(in, out, name)
AnalysisPoint(name) = AnalysisPoint(; name)

Base.show(io::IO, ap::AnalysisPoint) = show(io, MIME"text/plain"(), ap)
function Base.show(io::IO, ::MIME"text/plain", ap::AnalysisPoint)
    if get(io, :compact, false)
        print(io, "AnalysisPoint($(ap.in.u), $(ap.out.u); name=$(ap.name))")
    else
        print(io, "AnalysisPoint(")
        printstyled(io, ap.name, color = :cyan)
        if ap.in !== nothing && ap.out !== nothing
            print(io, " from ")
            printstyled(io, ap.in.u, color = :green)
            print(io, " to ")
            printstyled(io, ap.out.u, color = :blue)
        end
        print(io, ")")
    end
end

"""
    connect(in, ap::AnalysisPoint, out)
    connect(in, ap_name::Symbol, out)

Connect `in` and `out` with an [`AnalysisPoint`](@ref) inbetween.
The incoming connection `in` is expected to be of type [`RealOutput`](@ref), and vice versa.

# Arguments:
- `in`: A connector of type [`RealOutput`](@ref)
- `out`: A connector of type [`RealInput`](@ref)
- `ap`: An explicitly created [`AnalysisPoint`](@ref)
- `ap_name`: If a name is given, an [`AnalysisPoint`](@ref) with the given name will be created automatically.
"""
function ModelingToolkit.connect(in, ap::AnalysisPoint, out)
    ap.in = in
    ap.out = out
    return 0 ~ ap
end

function ModelingToolkit.connect(in, ap_name::Symbol, out)
    return 0 ~ AnalysisPoint(in, out, ap_name)
end

function ModelingToolkit.vars(ap::AnalysisPoint; op = Differential)
    vars(connect(ap.in, ap.out); op)
end

"""
    find_analysis_point(sys, name::Symbol)

Find and return the analysis point in `sys` with the specified `name`. If no matching [`AnalysisPoint`](@ref) is found, `nothing` is returned.
"""
function find_analysis_point(sys, name)
    sys = ModelingToolkit.flatten(sys)
    eqs = equations(sys)
    for eq in eqs
        eq.rhs isa AnalysisPoint && eq.rhs.name == name && (return eq.rhs)
    end
    nothing
end

"""
    expand_analysis_points(sys)

Replace analysis points with the identity connection connect(ap.in, ap.out). This is called before a system containing analysis points is simulated, in which case analysis points have no effect.
"""
function expand_analysis_points(sys)
    sys = ModelingToolkit.flatten(sys) # TODO: this does not namespace variables in connect statements properly https://github.com/SciML/ModelingToolkit.jl/issues/1826
    new_eqs = map(get_eqs(sys)) do eq
        eq.rhs isa AnalysisPoint || (return eq)
        ap = eq.rhs
        connect(ap.in, ap.out)
    end
    @set! sys.eqs = new_eqs
    sys
end

function ModelingToolkit.namespace_expr(ap::AnalysisPoint, sys, n = nameof(sys)) where {T}
    in = ModelingToolkit.renamespace(n, ap.in)
    out = ModelingToolkit.renamespace(n, ap.out)
    name = Symbol(n, :_, ap.name)
    AnalysisPoint(in, out, name)
end

function Base.:(==)(ap1::AnalysisPoint, ap2::AnalysisPoint)
    return ap1.in == ap2.in && ap1.out == ap2.out # Name doesn't really matter if inputs and outputs are the same
end

"""
    get_sensitivity(sys, ap::AnalysisPoint; kwargs)
    get_sensitivity(sys, ap_name::Symbol; kwargs)

Compute the sensitivity function in analysis point `ap`. The sensitivity function is obtained by introducing an infinitesimal perturbation `d` at the input of `ap`, linearizing the system and computing the transfer function between `d` and the output of `ap`.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref).
"""
function get_sensitivity(sys, ap::AnalysisPoint; kwargs...)
    sys = ModelingToolkit.flatten(sys) # To get namespacing right
    t = get_iv(sys)
    @variables d(t) = 0 # Perturbation serving as input to sensivity transfer function
    found = false
    new_eqs = map(equations(sys)) do eq
        eq.rhs == ap || (return eq)
        found = true
        ap.out.u ~ ap.in.u + d # This assumes that the connector has an internal vaiable named u
    end
    found || error("Did not find analysis point $ap")
    @set! sys.eqs = new_eqs
    @set! sys.states = [states(sys); d]
    @set! sys.defaults = merge(ModelingToolkit.defaults(sys), Dict(d => 0))
    sys = expand_analysis_points(sys) # Any remaining analysis points are removed by this
    ModelingToolkit.linearize(sys, [d], [ap.out.u]; kwargs...)
end

"""
    get_comp_sensitivity(sys, ap::AnalysisPoint; kwargs)
    get_comp_sensitivity(sys, ap_name::Symbol; kwargs)

Compute the complementary sensitivity function in analysis point `ap`. The complementary sensitivity function is obtained by introducing an infinitesimal perturbation `d` at the output of `ap`, linearizing the system and computing the transfer function between `d` and the input of `ap`.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_looptransfer`](@ref).
"""
function get_comp_sensitivity(sys, ap::AnalysisPoint; kwargs...)
    sys = ModelingToolkit.flatten(sys) # To get namespacing right
    t = get_iv(sys)
    @variables d(t) = 0 # Perturbation serving as input to sensivity transfer function
    found = false
    new_eqs = map(equations(sys)) do eq
        eq.rhs == ap || (return eq)
        found = true
        ap.out.u + d ~ ap.in.u # This assumes that the connector has an internal vaiable named u
    end
    found || error("Did not find analysis point $ap")
    @set! sys.eqs = new_eqs
    @set! sys.states = [states(sys); d]
    @set! sys.defaults = merge(ModelingToolkit.defaults(sys), Dict(d => 0))

    sys = expand_analysis_points(sys) # Any remaining analysis points are removed by this
    ModelingToolkit.linearize(sys, [d], [ap.in.u]; kwargs...)
end

"""
    get_looptransfer(sys, ap::AnalysisPoint; kwargs)
    get_looptransfer(sys, ap_name::Symbol; kwargs)

Compute the (linearized) loop-transfer function in analysis point `ap`, from `ap.out` to `ap.in`.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`open_loop`](@ref).
"""
function get_looptransfer(sys, ap::AnalysisPoint; kwargs...)
    sys = ModelingToolkit.flatten(sys) # To get namespacing right
    t = get_iv(sys)
    found = false
    new_eqs = map(equations(sys)) do eq
        eq.rhs == ap || (return eq)
        found = true
        0 ~ 0 # This assumes that the connector has an internal vaiable named u
    end
    found || error("Did not find analysis point $ap")
    @set! sys.eqs = new_eqs
    sys = expand_analysis_points(sys) # Any remaining analysis points are removed by this
    ModelingToolkit.linearize(sys, [ap.out.u], [ap.in.u]; kwargs...)
end

"""
    open_sys = open_loop(sys, ap::AnalysisPoint; kwargs)
    open_sys = open_loop(sys, ap_name::Symbol; kwargs)

Open the loop at analysis point `ap` by breaking the connection through `ap`.

`open_sys` will have `u ~ ap.out` as input and `y ~ ap.in` as output.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref).
"""
function open_loop(sys, ap::AnalysisPoint; kwargs...)
    sys = ModelingToolkit.flatten(sys) # To get namespacing right
    t = get_iv(sys)
    @variables u(t)=0 [input = true]
    @variables y(t)=0 [output = true]
    found = false
    new_eqs = map(equations(sys)) do eq
        eq.rhs == ap || (return [eq])
        found = true
        [ap.out.u ~ u
         ap.in.u ~ y]
    end
    found || error("Did not find analysis point $ap")
    new_eqs = reduce(vcat, new_eqs)
    @set! sys.eqs = new_eqs
    @set! sys.states = [states(sys); u; y]
    @set! sys.defaults = merge(ModelingToolkit.defaults(sys), Dict(u => 0, y => 0))
    sys
end

"""
    ModelingToolkit.linearize(sys, input::AnalysisPoint, output::AnalysisPoint)

Linearize a system between two analysis points.
All parts of the model that do not appear between `input` and `output` will be neglected.
"""
function ModelingToolkit.linearize(sys, input::AnalysisPoint, output::AnalysisPoint;
                                   kwargs...)
    sys = ModelingToolkit.flatten(sys) # To get namespacing right
    t = get_iv(sys)
    @variables u(t)=0 [input = true]
    @variables y(t)=0 [output = true]
    new_eqs = map(equations(sys)) do eq
        if eq.rhs == input
            [input.out.u ~ u]
            #input.in.u ~ 0] # We only need to ground one of the ends, hence not including this equation
        elseif eq.rhs == output
            [output.in.u ~ y
             output.out.u ~ 0]
        else
            return [eq]
        end
    end
    new_eqs = reduce(vcat, new_eqs)
    @set! sys.eqs = new_eqs
    @set! sys.states = [states(sys); u; y]
    @set! sys.defaults = merge(ModelingToolkit.defaults(sys), Dict(u => 0, y => 0))
    sys = expand_analysis_points(sys)
    ModelingToolkit.linearize(sys, [u], [y]; kwargs...)
end

# Add a method to get_sensitivity that accepts the name of an AnalysisPoint 
for f in [:get_sensitivity, :get_comp_sensitivity, :get_looptransfer, :open_loop]
    @eval function $f(sys, ap_name::Symbol, args...; kwargs...)
        ap = find_analysis_point(sys, ap_name)
        ap === nothing && error("Failed to find an analysis point named $ap_name")
        $f(sys, ap, args...; kwargs...)
    end
end

function ModelingToolkit.linearize(sys, input_name::Symbol, output_name::Symbol;
                                   kwargs...)
    input = find_analysis_point(sys, input_name)
    input === nothing && error("Failed to find an analysis point named $input_name")
    output = find_analysis_point(sys, output_name)
    output === nothing && error("Failed to find an analysis point named $output_name")
    ModelingToolkit.linearize(sys, input, output; kwargs...)
end
