# Running Models with Discrete Data

There are 3 ways to include data as part of a model.

 1. using `ModelingToolkitStandardLibrary.Blocks.TimeVaryingFunction`
 2. using a custom component with external data
 3. using `ModelingToolkitStandardLibrary.Blocks.SampledData`

This tutorial demonstrate each case and explain the pros and cons of each.

## `TimeVaryingFunction` Component

The `ModelingToolkitStandardLibrary.Blocks.TimeVaryingFunction` component is easy to use and is performant.  However the data is locked to the `ODESystem` and can only be changed by building a new `ODESystem`.  Therefore, running a batch of data would not be efficient.  Below is an example of how to use the `TimeVaryingFunction` with `DataInterpolations` to build the function from sampled discrete data.

```julia
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using DataInterpolations
using OrdinaryDiffEq

@parameters t
D = Differential(t)

function System(f; name)
    src = TimeVaryingFunction(f)

    vars = @variables f(t)=0 x(t)=0 dx(t)=0 ddx(t)=0
    pars = @parameters m=10 k=1000 d=1

    eqs = [f ~ src.output.u
        ddx * 10 ~ k * x + d * dx + f
        D(x) ~ dx
        D(dx) ~ ddx]

    ODESystem(eqs, t, vars, pars; systems = [src], name)
end

dt = 4e-4
time = 0:dt:0.1
data = sin.(2 * pi * time * 100) # example data

f = LinearInterpolation(data, time)

@named system = System(f)
sys = structural_simplify(system)
prob = ODEProblem(sys, [], (0, time[end]))
sol = solve(prob, ImplicitEuler())
```

If we want to run a new data set, this requires building a new `LinearInterpolation` and `ODESystem` followed by running `structural_simplify`, all of which takes time.  Therefore, to run several pieces of data it's better to re-use an `ODESystem`.  The next couple methods will demonstrate how to do this.

## Custom Component with External Data

The below code shows how to include data using a `Ref` and registered `get_sampled_data` function.  This example uses a very basic function which requires non-adaptive solving and sampled data.  As can be seen, the data can easily be set and changed before solving.

```julia
const rdata = Ref{Vector{Float64}}()

# Data Sets
data1 = sin.(2 * pi * time * 100)
data2 = cos.(2 * pi * time * 50)

function get_sampled_data(t)
    i = floor(Int, t / dt) + 1
    x = rdata[][i]

    return x
end

Symbolics.@register_symbolic get_sampled_data(t)

function System(; name)
    vars = @variables f(t)=0 x(t)=0 dx(t)=0 ddx(t)=0
    pars = @parameters m=10 k=1000 d=1

    eqs = [f ~ get_sampled_data(t)
        ddx * 10 ~ k * x + d * dx + f
        D(x) ~ dx
        D(dx) ~ ddx]

    ODESystem(eqs, t, vars, pars; name)
end

@named system = System()
sys = structural_simplify(system)
prob = ODEProblem(sys, [], (0, time[end]))

rdata[] = data1
sol1 = solve(prob, ImplicitEuler(); dt, adaptive = false)
ddx1 = sol1[sys.ddx]

rdata[] = data2
sol2 = solve(prob, ImplicitEuler(); dt, adaptive = false)
ddx2 = sol2[sys.ddx]
```

The drawback of this method is that the solution observables can be linked to the data `Ref`, which means that if the data changes then the observables are no longer valid.  In this case `ddx` is an observable that is derived directly from the data.  Therefore, `sol1[sys.ddx]` is no longer correct after the data is changed for `sol2`.

```julia
# the following test will fail
@test all(ddx1 .== sol1[sys.ddx]) #returns false
```

Additional code could be added to resolve this issue, for example by using a `Ref{Dict}` that could link a parameter of the model to the data source.  This would also be necessary for parallel processing.

## `SampledData` Component

To resolve the issues presented above, the `ModelingToolkitStandardLibrary.Blocks.SampledData` component can be used which allows for a resusable `ODESystem` and self contained data which ensures a solution which remains valid for it's lifetime.  Now it's possible to also parallelize the call to `solve()`.

```julia
function System(; name)
    src = SampledData(Float64)

    vars = @variables f(t)=0 x(t)=0 dx(t)=0 ddx(t)=0
    pars = @parameters m=10 k=1000 d=1

    eqs = [f ~ src.output.u
        ddx * 10 ~ k * x + d * dx + f
        D(x) ~ dx
        D(dx) ~ ddx]

    ODESystem(eqs, t, vars, pars; systems = [src], name)
end

@named system = System()
sys = structural_simplify(system)
s = complete(system)
prob = ODEProblem(sys, [], (0, time[end]); tofloat = false)
defs = ModelingToolkit.defaults(sys)

function get_prob(data)
    defs[s.src.buffer] = Parameter(data, dt)
    # ensure p is a uniform type of Vector{Parameter{Float64}} (converting from Vector{Any})
    p = Parameter.(ModelingToolkit.varmap_to_vars(defs, parameters(sys); tofloat = false))
    remake(prob; p)
end

prob1 = get_prob(data1)
prob2 = get_prob(data2)

sol1 = Ref{ODESolution}()
sol2 = Ref{ODESolution}()
@sync begin
    @async sol1[] = solve(prob1, ImplicitEuler())
    @async sol2[] = solve(prob2, ImplicitEuler())
end
```

Note, in the above example, we can build the system with an empty `SampledData` component, only setting the expected data type: `@named src = SampledData(Float64)`.  It's also possible to initialize the component with real sampled data: `@named src = SampledData(data, dt)`.  Additionally note that before running an `ODEProblem` using the `SampledData` component, one must be careful about the parameter vector Type.  The `SampledData` component contains a `buffer` parameter of type `Parameter`, therefore we must generate the problem using `tofloat=false`.  This will initially give a parameter vector of type `Vector{Any}` with a mix of numbers and `Parameter` type.  We can convert the vector to a uniform `Parameter` type by running `p = Parameter.(p)`.  This will wrap all the single values in a `Parameter` which will be mathematically equivalent to a `Number`.
