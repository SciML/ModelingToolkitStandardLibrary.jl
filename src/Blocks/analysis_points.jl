using ModelingToolkit: get_eqs, vars, @set!, get_iv

Base.@kwdef mutable struct AnalysisPoint
    in = nothing
    out = nothing
    name::Symbol = :nothing
end

Base.nameof(ap::AnalysisPoint) = ap.name

"""
    AnalysisPoint(in, out, name::Symbol)
    AnalysisPoint(in, out; name::Symbol)
    AnalysisPoint(name::Symbol)

Create an AnalysisPoint for linear analysis. Analysis points can also be created automatically by calling
```
connect(in, :ap_name, out)
```

!!! danger "Experimental"
    The analysis-point interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning.

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
    if ap.in === nothing
        print(io, "0")
        return
    end
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

_isinput(x) = x isa ODESystem && endswith(string(nameof(x)), "input")
_isoutput(x) = x isa ODESystem && endswith(string(nameof(x)), "output")
function ap_warning(n)
    @warn "The $(n == 1 ? "first" : "third") argument to a connection with an analysis point was a $(n == 1 ? "RealInput" : "RealOutput"). This is supported in order to handle inverse models, but may not be what you intended. If you are building a forward model (causal), you may want to swap the first and the third arguments to connect. Learn more about the causality of analysis points in the docstring for AnalysisPoint. Silence this message by connect(out, :name, in; verbose = false)"
end

"""
    connect(output_connector, ap_name::Symbol, input_connector; verbose = true)
    connect(output_connector, ap::AnalysisPoint, input_connector; verbose = true)

Connect `output_connector` and `input_connector` with an [`AnalysisPoint`](@ref) inbetween.
The incoming connection `output_connector` is expected to be of type [`RealOutput`](@ref), and vice versa.
NOTE: The connection is assumed to be *causal*, meaning that
```
connect(C.output, :plant_input, P.input)
```
is correct, whereas
```
connect(P.input, :plant_input, C.output)
```
typically is not (unless the model is an inverse model).

# Arguments:
- `output_connector`: A connector of type [`RealOutput`](@ref)
- `input_connector`: A connector of type [`RealInput`](@ref)
- `ap`: An explicitly created [`AnalysisPoint`](@ref)
- `ap_name`: If a name is given, an [`AnalysisPoint`](@ref) with the given name will be created automatically.
- `verbose`: Causes a warning to be displayed if an input is connected to an output (reverse causality). Silence this warning if you are analysing an inverse model.
"""
function ModelingToolkit.connect(in, ap::AnalysisPoint, out; verbose = true)
    verbose && _isinput(in) && ap_warning(1)
    verbose && _isoutput(out) && ap_warning(2)
    ap.in = in
    ap.out = out
    return AnalysisPoint() ~ ap
end

function ModelingToolkit.connect(in, ap_name::Symbol, out; verbose = true)
    verbose && _isinput(in) && ap_warning(1)
    verbose && _isoutput(out) && ap_warning(2)
    ap = AnalysisPoint(in, out, ap_name)
    return AnalysisPoint() ~ ap
end

ModelingToolkit.get_systems(ap::AnalysisPoint) = (ap.in, ap.out)

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

function Base.:(==)(ap1::AnalysisPoint, ap2::AnalysisPoint)
    return ap1.in == ap2.in && ap1.out == ap2.out # Name doesn't really matter if inputs and outputs are the same
end

"""
    get_sensitivity(sys, ap::AnalysisPoint; kwargs)
    get_sensitivity(sys, ap_name::Symbol; kwargs)

Compute the sensitivity function in analysis point `ap`. The sensitivity function is obtained by introducing an infinitesimal perturbation `d` at the input of `ap`, linearizing the system and computing the transfer function between `d` and the output of `ap`.

!!! danger "Experimental"
    The analysis-point interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref).
"""
function get_sensitivity(sys, ap_name::Symbol; kwargs...)
    find = function (x, ns)
        x isa AnalysisPoint || return false
        if ns === nothing
            nameof(x) === ap_name
        else
            Symbol(ns, :_, nameof(x)) === ap_name
        end
    end
    t = get_iv(sys)
    @variables d(t) = 0 # Perturbation serving as input to sensitivity transfer function
    namespace = Ref{Union{Nothing, Symbol}}(nothing)
    apr = Ref{Union{Nothing, AnalysisPoint}}(nothing)
    replace = let d = d, namespace = namespace, apr = apr
        function (ap, ns)
            namespace[] = ns # Save the namespace to make it available for renamespace below
            apr[] = ap
            (ap.out.u ~ ap.in.u + d), d
        end
    end
    sys = expand_connections(sys, find, replace)
    (ap = apr[]) === nothing && error("Did not find analysis point $ap")
    u = ap.out.u
    if (ns = namespace[]) !== nothing
        d = ModelingToolkit.renamespace(ns, d)
        u = ModelingToolkit.renamespace(ns, u)
    end
    ModelingToolkit.linearize(sys, [d], [u]; kwargs...)
end

"""
    get_comp_sensitivity(sys, ap::AnalysisPoint; kwargs)
    get_comp_sensitivity(sys, ap_name::Symbol; kwargs)

Compute the complementary sensitivity function in analysis point `ap`. The complementary sensitivity function is obtained by introducing an infinitesimal perturbation `d` at the output of `ap`, linearizing the system and computing the transfer function between `d` and the input of `ap`.

!!! danger "Experimental"
    The analysis-point interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_looptransfer`](@ref).
"""
function get_comp_sensitivity(sys, ap_name::Symbol; kwargs...)
    find = function (x, ns)
        x isa AnalysisPoint || return false
        if ns === nothing
            nameof(x) === ap_name
        else
            Symbol(ns, :_, nameof(x)) === ap_name
        end
    end
    t = get_iv(sys)
    @variables d(t) = 0 # Perturbation serving as input to sensitivity transfer function
    namespace = Ref{Union{Nothing, Symbol}}(nothing)
    apr = Ref{Union{Nothing, AnalysisPoint}}(nothing)
    replace = let d = d, namespace = namespace, apr = apr
        function (ap, ns)
            namespace[] = ns # Save the namespace to make it available for renamespace below
            apr[] = ap
            (ap.out.u + d ~ ap.in.u), d
        end
    end
    sys = expand_connections(sys, find, replace)
    (ap = apr[]) === nothing && error("Did not find analysis point $ap")
    u = ap.in.u
    if (ns = namespace[]) !== nothing
        d = ModelingToolkit.renamespace(ns, d)
        u = ModelingToolkit.renamespace(ns, u)
    end
    ModelingToolkit.linearize(sys, [d], [u]; kwargs...)
end

"""
    get_looptransfer(sys, ap::AnalysisPoint; kwargs)
    get_looptransfer(sys, ap_name::Symbol; kwargs)

Compute the (linearized) loop-transfer function in analysis point `ap`, from `ap.out` to `ap.in`.

!!! danger "Experimental"
    The analysis-point interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`open_loop`](@ref).
"""
function get_looptransfer(sys, ap_name::Symbol; kwargs...)
    find = function (x, ns)
        x isa AnalysisPoint || return false
        if ns === nothing
            nameof(x) === ap_name
        else
            Symbol(ns, :_, nameof(x)) === ap_name
        end
    end
    t = get_iv(sys)
    namespace = Ref{Union{Nothing, Symbol}}(nothing)
    apr = Ref{Union{Nothing, AnalysisPoint}}(nothing)
    replace = let namespace = namespace, apr = apr
        function (ap, ns)
            namespace[] = ns # Save the namespace to make it available for renamespace below
            apr[] = ap
            (0 ~ 0), nothing
        end
    end
    sys = expand_connections(sys, find, replace)
    (ap = apr[]) === nothing && error("Did not find analysis point $ap")
    u = ap.out.u
    y = ap.in.u
    if (ns = namespace[]) !== nothing
        y = ModelingToolkit.renamespace(ns, y)
        u = ModelingToolkit.renamespace(ns, u)
    end
    ModelingToolkit.linearize(sys, [u], [y]; kwargs...)
end

"""
    open_sys = open_loop(sys, ap::AnalysisPoint; kwargs)
    open_sys = open_loop(sys, ap_name::Symbol; kwargs)

Open the loop at analysis point `ap` by breaking the connection through `ap`.

`open_sys` will have `u ~ ap.out` as input and `y ~ ap.in` as output.

!!! danger "Experimental"
    The analysis-point interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning.

# Arguments:
- `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`get_looptransfer`](@ref).
"""
function open_loop(sys, ap_name::Symbol; kwargs...)
    find = function (x, ns)
        x isa AnalysisPoint || return false
        if ns === nothing
            nameof(x) === ap_name
        else
            Symbol(ns, :_, nameof(x)) === ap_name
        end
    end
    t = get_iv(sys)
    @variables u(t)=0 [input = true]
    @variables y(t)=0 [output = true]
    namespace = Ref{Union{Nothing, Symbol}}(nothing)
    apr = Ref{Union{Nothing, AnalysisPoint}}(nothing)
    replace = let u = u, y = y, namespace = namespace, apr = apr
        function (ap, ns)
            namespace[] = ns # Save the namespace to make it available for renamespace below
            apr[] = ap
            [ap.out.u ~ u, ap.in.u ~ y], [u, y]
        end
    end
    if (ns = namespace[]) !== nothing
        y = ModelingToolkit.renamespace(ns, y)
        u = ModelingToolkit.renamespace(ns, u)
    end
    sys = expand_connections(sys, find, replace)
    (ap = apr[]) === nothing && error("Did not find analysis point $ap")
    sys
end

"""
    ModelingToolkit.linearize(sys, input_name::Symbol, output_name::Symbol)

Linearize a system between two analysis points. To get a loop-transfer function, see [`get_looptransfer`](@ref)
"""
function ModelingToolkit.linearize(sys, input_name::Symbol, output_name::Symbol;
                                   kwargs...)
    find = function (x, ns)
        x isa AnalysisPoint || return false
        if ns === nothing
            nameof(x) ∈ (input_name, output_name)
        else
            Symbol(ns, :_, nameof(x)) ∈ (input_name, output_name)
        end
    end
    t = get_iv(sys)
    @variables u(t)=0 [input = true]
    @variables y(t)=0 [output = true]
    namespace = Ref{Union{Nothing, Symbol}}(nothing)
    apr = Ref{Union{Nothing, AnalysisPoint}}(nothing)
    replace = let u = u, y = y, namespace = namespace, apr = apr
        function (ap, ns)
            namespace[] = ns # Save the namespace to make it available for renamespace below
            apr[] = ap
            if nameof(ap) === input_name
                [ap.out.u ~ ap.in.u + u], u
                #input.in.u ~ 0] # We only need to ground one of the ends, hence not including this equation
            elseif nameof(ap) === output_name
                [ap.in.u ~ y
                 ap.out.u ~ ap.in.u], y
            else
                error("This should never happen")
            end
        end
    end
    sys = expand_connections(sys, find, replace)
    (ap = apr[]) === nothing && error("Did not find analysis point $ap")
    if (ns = namespace[]) !== nothing
        y = ModelingToolkit.renamespace(ns, y)
        u = ModelingToolkit.renamespace(ns, u)
    end
    ModelingToolkit.linearize(sys, [u], [y]; kwargs...)
end

# Add a method to get_sensitivity that accepts the name of an AnalysisPoint
for f in [:get_sensitivity, :get_comp_sensitivity, :get_looptransfer, :open_loop]
    @eval function $f(sys, ap::AnalysisPoint, args...; kwargs...)
        $f(sys, nameof(ap), args...; kwargs...)
    end
end

function ModelingToolkit.linearize(sys, input::AnalysisPoint, output::AnalysisPoint;
                                   kwargs...)
    ModelingToolkit.linearize(sys, nameof(input), nameof(output); kwargs...)
end
