using ModelingToolkit: get_eqs, vars, @set!, get_iv

Base.@kwdef mutable struct AnalysisPoint
    in = nothing
    out = nothing
    name::Symbol = :nothing
end
Base.broadcastable(x::AnalysisPoint) = Ref(x)
if Base.isdefined(ModelingToolkit, :isconnection)
    ModelingToolkit.isconnection(::AnalysisPoint) = true
end

Base.nameof(ap::AnalysisPoint) = ap.name
function Base.hash(ap::AnalysisPoint, seed::UInt)
    h1 = hash(ap.in, seed)
    h2 = hash(ap.out, h1)
    h3 = hash(ap.name, h2)
    h3 ⊻ (0xd29cdc51aa6562d4 % UInt)
end

function ModelingToolkit.get_unit(ap::AnalysisPoint)
    ModelingToolkit.unitless
end

function ap_var(sys)
    if hasproperty(sys, :u)
        # collect to turn symbolic arrays into arrays of symbols
        return length(sys.u) == 1 ? sys.u : collect(sys.u)
    end
    x = unknowns(sys)
    length(x) == 1 && return x[1]
    error("Could not determine the analysis-point variable in system $(nameof(sys)). To use an analysis point, apply it to a connection between two causal blocks containing connectors of type `RealInput/RealOutput` from ModelingToolkitStandardLibrary.Blocks.")
end

"""
    find_analysis_points(sys)

Return a list of all analysis points in `sys`. If none are found, the list is empty.
"""
function find_analysis_points(sys)
    sys = ModelingToolkit.flatten(sys)
    eqs = equations(sys)
    aps = []
    for eq in eqs
        if eq.rhs isa AnalysisPoint
            push!(aps, eq.rhs)
        end
    end
    aps
end

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
@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(; k = -1)
t = ModelingToolkit.get_iv(P)
eqs = [connect(P.output, C.input)
       connect(C.output, :plant_input, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :feedback_system)

matrices_S, _ = get_sensitivity(sys, :plant_input) # Compute the matrices of a state-space representation of the (input) sensitivity function.
matrices_T, _ = get_comp_sensitivity(sys, :plant_input)
```

Continued linear analysis and design can be performed using ControlSystemsBase.jl.
Create `ControlSystemsBase.StateSpace` objects using

```julia
using ControlSystemsBase, Plots
S = ss(matrices_S...)
T = ss(matrices_T...)
bodeplot([S, T], lab = ["S" "T"])
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
        print(io, "AnalysisPoint($(ap_var(ap.in)), $(ap_var(ap.out)); name=$(ap.name))")
    else
        print(io, "AnalysisPoint(")
        printstyled(io, ap.name, color = :cyan)
        if ap.in !== nothing && ap.out !== nothing
            print(io, " from ")
            printstyled(io, ap_var(ap.in), color = :green)
            print(io, " to ")
            printstyled(io, ap_var(ap.out), color = :blue)
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

*PLEASE NOTE*: The connection is assumed to be *causal*, meaning that

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

function namespaced_ap_match(ap_name, loop_openings)
    (x, ns) -> namespaced_ap_match(x, ns, ap_name, loop_openings)
end # 🍛

"""
    namespaced_ap_match(x, ns, ap_names, loop_openings)

Returns true if `x` is an AnalysisPoint and matches either any of `ap_names` or any of the loop openings if namedspaced with `ns`.
"""
function namespaced_ap_match(x, ns, ap_names0, loop_openings)
    x isa AnalysisPoint || return false
    ap_names = Set(ap_names0 isa Symbol ? [ap_names0] : ap_names0)
    if ns === nothing
        nameof(x) ∈ ap_names || (loop_openings !== nothing && nameof(x) ∈ loop_openings)
    else
        xx = Symbol(ns, :_, nameof(x))
        xx ∈ ap_names || (loop_openings !== nothing && xx ∈ loop_openings)
    end
end

function get_perturbation_var(x::Num, prefix = "d")
    @variables d(t) = 0
    @set! d.val.f.name = Symbol("$(prefix)_$(x)")
    d
end
function get_perturbation_var(x, args...)
    get_perturbation_var.(x, args...)
end

function _check_and_sort!(ap_names, aps, namespaces, multiplicities)
    ap_names isa Symbol && (ap_names = [ap_names])
    happy = true
    permutation = Int[]
    for apn in ap_names
        ind = findfirst(eachindex(aps)) do i
            x = aps[i]
            ns = namespaces[i]
            ns === nothing ? nameof(x) == apn : Symbol(ns, :_, nameof(x)) == apn
        end
        if ind === nothing
            @error "Could not find analysis point $apn"
            happy = false
        end
        push!(permutation, ind)
    end
    happy || error("Failed to find all analysis points. I found these: $(nameof.(aps))")
    aps .= aps[permutation]
    namespaces .= namespaces[permutation]
    # extend permutation to account for aps that introduce several vars (array-valued connect statements)
    # This requires two passes, the first adjusts the permutation indices, the second inserts the new permutation indices
    for i in eachindex(permutation)
        pᵢ = permutation[i]
        mi = multiplicities[i]
        if mi > 1
            # all permutations that have a permutation index higher than pᵢ need to be incremented with mi-1
            # since we are inserting mi-1 new permutation indices
            for j in eachindex(permutation)
                pⱼ = permutation[j]
                pⱼ > pᵢ && (permutation[j] = pⱼ + mi - 1)
            end
        end
    end
    new_perm = Int[]
    for i in eachindex(permutation)
        pᵢ = permutation[i]
        mi = multiplicities[i]
        if mi > 1
            new_pi = (1:mi) .+ (pᵢ - 1) # these are the new permutations arising from apᵢ
            append!(new_perm, new_pi)
        else
            push!(new_perm, pᵢ)
        end
    end
    @assert length(new_perm) == sum(multiplicities)
    @assert sort(new_perm) == 1:length(new_perm)
    new_perm
end

const SymOrVec = Union{Symbol, Vector{Symbol}}

function get_sensitivity_function(
        sys, ap_name::SymOrVec; loop_openings = nothing, system_modifier = identity,
        kwargs...)
    find = namespaced_ap_match(ap_name, loop_openings)
    t = get_iv(sys)
    aps = []
    u = []
    d = []
    multiplicities = Int[]
    namespaces = []
    replace = function (ap, ns)
        if namespaced_ap_match(ap, ns, ap_name, nothing)
            di = get_perturbation_var(ap_var(ap.in))
            push!(aps, ap)
            push!(multiplicities, length(di)) # one ap may yield several new vars
            push!(namespaces, ns)
            append!(d, di)
            append!(u, ap_var(ap.out))
            (ap_var(ap.out) .~ ap_var(ap.in) + di), di
        else # loop opening
            [ap_var(ap.out) .~ 0;], []
        end
    end
    sys = expand_connections(sys, find, replace)
    permutation = _check_and_sort!(ap_name, aps, namespaces, multiplicities)
    dn = ModelingToolkit.renamespace.(namespaces, d[permutation])
    un = ModelingToolkit.renamespace.(namespaces, u[permutation])
    sys = system_modifier(sys)
    ModelingToolkit.linearization_function(sys, dn, un; kwargs...)
end

function get_comp_sensitivity_function(
        sys, ap_name::SymOrVec; loop_openings = nothing, system_modifier = identity,
        kwargs...)
    find = namespaced_ap_match(ap_name, loop_openings)
    t = get_iv(sys)
    aps = []
    u = []
    d = []
    multiplicities = Int[]
    namespaces = []
    replace = function (ap, ns)
        if namespaced_ap_match(ap, ns, ap_name, nothing)
            di = get_perturbation_var(ap_var(ap.in))
            push!(aps, ap)
            push!(multiplicities, length(di)) # one ap may yield several new vars
            push!(namespaces, ns)
            append!(d, di)
            append!(u, ap_var(ap.in))
            (ap_var(ap.out) + di .~ ap_var(ap.in)), di
        else # loop opening
            [ap_var(ap.out) .~ 0;], []
        end
    end
    sys = expand_connections(sys, find, replace)
    permutation = _check_and_sort!(ap_name, aps, namespaces, multiplicities)
    dn = ModelingToolkit.renamespace.(namespaces, d[permutation])
    un = ModelingToolkit.renamespace.(namespaces, u[permutation])
    sys = system_modifier(sys)
    ModelingToolkit.linearization_function(sys, dn, un; kwargs...)
end

function get_looptransfer_function(
        sys, ap_name::SymOrVec; loop_openings = nothing, system_modifier = identity,
        kwargs...)
    find = namespaced_ap_match(ap_name, loop_openings)
    t = get_iv(sys)
    aps = []
    multiplicities = Int[]
    namespaces = []
    replace = function (ap, ns)
        if namespaced_ap_match(ap, ns, ap_name, nothing)
            push!(aps, ap)
            push!(multiplicities, length(ap_var(ap.in)))
            push!(namespaces, ns)
            (0 ~ 0), nothing
        else # loop opening
            [ap_var(ap.out) .~ 0;], []
        end
    end
    sys = expand_connections(sys, find, replace)
    permutation = _check_and_sort!(ap_name, aps, namespaces, multiplicities)
    u = reduce(vcat, ap_var(ap.out) for ap in aps)
    y = reduce(vcat, ap_var(ap.in) for ap in aps)
    yn = ModelingToolkit.renamespace.(namespaces, y)# permutation applied in _check_and_sort
    un = ModelingToolkit.renamespace.(namespaces, u)
    sys = system_modifier(sys)
    ModelingToolkit.linearization_function(sys, un, yn; kwargs...)
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
function open_loop(
        sys, ap_name::Symbol; ground_input = false, system_modifier = identity, kwargs...)
    find = namespaced_ap_match(ap_name, nothing)
    t = get_iv(sys)
    @variables u(t)=0 [input = true]
    @variables y(t)=0 [output = true]
    namespace = Ref{Union{Nothing, Symbol}}(nothing)
    apr = Ref{Union{Nothing, AnalysisPoint}}(nothing)
    replace = let u = u, y = y, namespace = namespace, apr = apr
        function (ap, ns)
            namespace[] = ns # Save the namespace to make it available for renamespace below
            apr[] = ap
            if ground_input
                [ap_var(ap.out) ~ 0, ap_var(ap.in) ~ y], [y]
            else
                [ap_var(ap.out) ~ u, ap_var(ap.in) ~ y], [u, y]
            end
        end
    end
    sys = expand_connections(sys, find, replace)
    (ap = apr[]) === nothing && error("Did not find analysis point $ap_name")
    sys = system_modifier(sys)
end

function ModelingToolkit.linearization_function(sys::ModelingToolkit.AbstractSystem,
        input_name::SymOrVec, output_name;
        loop_openings = nothing,
        system_modifier = identity, # This is used to, e.g., apply JuliaSimCompiler.IRSystem after analysis-point handling
        kwargs...)
    t = get_iv(sys)
    @variables u(t)=0 [input = true]
    names = [input_name;]
    if output_name isa SymOrVec
        @variables y(t)=0 [output = true]
        names = [names; output_name]
    end
    find = namespaced_ap_match(names, loop_openings)

    u = []
    y = []
    namespace_u = []
    namespace_y = []
    aps_u = []
    aps_y = []
    multiplicities_u = Int[]
    multiplicities_y = Int[]

    replace = function (ap, ns)
        if namespaced_ap_match(ap, ns, input_name, nothing)
            push!(namespace_u, ns) # Save the namespace to make it available for renamespace below
            push!(aps_u, ap)
            ui = get_perturbation_var(ap_var(ap.out), "u")
            push!(multiplicities_u, length(ui)) # one ap may yield several new vars
            append!(u, ui)
            if loop_openings !== nothing && ap.name ∈ loop_openings
                # In this case, we break the existing connection.
                [ap_var(ap.out) .~ ui;], ui
            else
                [ap_var(ap.out) .~ ap_var(ap.in) + ui;], ui
            end
            #input.in.u ~ 0] # We only need to ground one of the ends, hence not including this equation
        elseif output_name isa SymOrVec && namespaced_ap_match(ap, ns, output_name, nothing)
            push!(namespace_y, ns) # Save the namespace to make it available for renamespace below
            push!(aps_y, ap)
            yi = get_perturbation_var(ap_var(ap.in), "y")
            push!(multiplicities_y, length(yi))
            append!(y, yi)
            if loop_openings !== nothing && ap.name ∈ loop_openings
                [ap_var(ap.in) .~ yi;
                 ap_var(ap.out) .~ 0], yi # In this case, we break the existing connection.
            else
                [ap_var(ap.in) .~ yi;
                 ap_var(ap.out) .~ ap_var(ap.in)], yi
            end
        else # loop opening
            [ap_var(ap.out) .~ 0;], []
        end
    end

    sys = expand_connections(sys, find, replace)

    permutation_u = _check_and_sort!(input_name, aps_u, namespace_u, multiplicities_u)
    un = ModelingToolkit.renamespace.(namespace_u, u)
    if output_name isa SymOrVec
        permutation_y = _check_and_sort!(output_name, aps_y, namespace_y, multiplicities_y)
        yn = ModelingToolkit.renamespace.(namespace_y, y) # permutation applied in _check_and_sort
    else
        yn = output_name
    end
    ModelingToolkit.linearization_function(system_modifier(sys), un, yn; kwargs...)
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

# Methods above are implemented in terms of linearization_function, the method below creates wrappers for linearize
for f in [:get_sensitivity, :get_comp_sensitivity, :get_looptransfer]
    @eval function $f(sys,
            ap,
            args...;
            loop_openings = nothing,
            op = Dict(),
            p = DiffEqBase.NullParameters(),
            system_modifier = identity,
            kwargs...)
        lin_fun, ssys = $(Symbol(string(f) * "_function"))(sys, ap, args...; op, p,
            loop_openings, system_modifier,
            kwargs...)
        ModelingToolkit.linearize(ssys, lin_fun; op, p, kwargs...), ssys
    end
end

"""
    ModelingToolkit.linearize(sys, input_name::Symbol, output_name; kwargs...)

Linearize a system between two analysis points. To get a loop-transfer function, see [`get_looptransfer`](@ref).

The output is allowed to be either an analysis-point name, or a vector of symbolic variables like the standard interface to `linearize`. The input must be an analysis-point name.
"""
function ModelingToolkit.linearize(sys, input_name::SymOrVec, output_name;
        loop_openings = nothing, system_modifier = identity, kwargs...)
    lin_fun, ssys = linearization_function(sys, input_name, output_name;
        loop_openings, system_modifier, kwargs...)
    ModelingToolkit.linearize(ssys, lin_fun; kwargs...), ssys
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
get_sensitivity

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
get_comp_sensitivity

"""
    get_looptransfer(sys, ap::AnalysisPoint; kwargs)
    get_looptransfer(sys, ap_name::Symbol; kwargs)

Compute the (linearized) loop-transfer function in analysis point `ap`, from `ap.out` to `ap.in`.

!!! info "Negative feedback"

    Feedback loops often use negative feedback, and the computed loop-transfer function will in this case have the negative feedback included. Standard analysis tools often assume a loop-transfer function without the negative gain built in, and the result of this function may thus need negation before use.


!!! danger "Experimental"

    The analysis-point interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning.

# Arguments:

  - `kwargs`: Are sent to `ModelingToolkit.linearize`

See also [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref), [`open_loop`](@ref).
"""
get_looptransfer
