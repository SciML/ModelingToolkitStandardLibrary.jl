# Cauer Low Pass Analog Filter

The example Cauer Filter is a low-pass-filter of the fifth order.
It is realized using an analog network.
The [`Electrical.Voltage`](@ref) source is the input voltage (whose value varies as defined by [`Blocks.Step`](@ref)), and the `resistor.p.v` is the filter output voltage.
The pulse response is calculated.

## Copy-Pastable Example

```@example cauer_low_pass_analog
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks: Step
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots

@variables t

@mtkmodel CauerLowPassAnalog begin
    @parameters begin
        L1 = 1.304,             [description="Inductance filter coefficient no.1"]
        L2 = 0.8586,            [description="Inductance filter coefficient no.2"]
        C1 = 1.072,             [description="Capacitance filter coefficient no.1"]
        C2 = 1/(1.704992^2*L1), [description="Capacitance filter coefficient no.2"]
        C3 = 1.682,             [description="Capacitance filter coefficient no.3"]
        C4 = 1/(1.179945^2*L2), [description="Capacitance filter coefficient no.4"]
        C5 = 0.7262,            [description="Capacitance filter coefficient no.5"]
    end
    @components begin
        step = Step(height=1, start_time=1, smooth=false)
        source = Voltage()
        ground = Ground()
        resistor1 = Resistor(R=1)
        resistor2 = Resistor(R=1)
        inductor1 = Inductor(L=L1)
        inductor2 = Inductor(L=L2)
        capacitor1 = Capacitor(C=C1)
        capacitor2 = Capacitor(C=C2)
        capacitor3 = Capacitor(C=C3)
        capacitor4 = Capacitor(C=C4)
        capacitor5 = Capacitor(C=C5)
    end
    @equations begin
        connect(step.output, source.V)
        connect(resistor1.n, capacitor1.p)
        connect(capacitor1.n, ground.g)
        connect(inductor1.p, capacitor2.p)
        connect(inductor1.p, capacitor1.p)
        connect(inductor1.n, capacitor2.n)
        connect(capacitor2.n, capacitor3.p)
        connect(capacitor2.n, capacitor4.p)
        connect(capacitor2.n, inductor2.p)
        connect(inductor2.n, capacitor4.n)
        connect(capacitor4.n, capacitor5.p)
        connect(capacitor4.n, resistor2.p)
        connect(capacitor1.n, capacitor3.n)
        connect(capacitor1.n, capacitor5.n)
        connect(resistor2.n, capacitor1.n)
        connect(resistor1.p, source.p)
        connect(source.n, ground.g)
    end
end

@mtkbuild model = CauerLowPassAnalog()

tspan = (0.0, 60.0)
prob = ODEProblem(model, ModelingToolkit.missing_variable_defaults(model), tspan)
sol = solve(prob, Rosenbrock23());
Plots.plot(sol; idxs=[model.source.p.v, model.resistor2.p.v])
Plots.savefig("cauer_low_pass_analog.png"); nothing # hide
```

![Cauer Low Pass Analog Filter plot](cauer_low_pass_analog.png)

## Explanation

### Setting up the Environment

Each component needed for this example is defined in the [Electrical Components](@ref "ModelingToolkitStandardLibrary: Electrical Components") module with the exception of [`Blocks.Step`](@ref).
These modules are loaded along with 

```julia
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks: Step
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots
```

### Defining Independent Variable

As usual, you must specify the independent variable time, `t`.
If prefered, you may alternatively run `using ModelingToolkitStandardLibrary.Blocks.t` to load [`Blocks.t`].

```julia
@variables t
```

### Defining the Model

This example uses the `@mtkmodel` macro to define the model.
This is the recommended method for defining models using `ModelingToolkit` and
provides several conveniences by way of reduced boilerplate when compared to past methods.
You can name the model and prepare the scaffold which you'll fill in below.

```julia
@mtkmodel CauerLowPassAnalog begin
    @parameters begin #= ... =# end
    @components begin #= ... =# end
    @equations begin #= ... =# end
end
```

#### Defining Model Parameters

There are a few interesting aspects to how the parameters are defined.

First, note that you are able to define parameters in terms of other paramters.
In this example, capacitance coefficients `C2` and `C4` are defined in terms
of inductance coefficients `L1` and `L2` respectively.

Second, note that you may optionally provide descriptions for each parameter.
The descriptive string is stored with the parameter and can be easily retrieve e.g., by `ModelingToolkit.getdescription(L1)`.
Storing and retrieving metadata can be useful as your team and models grow over time.
The descriptive string acts as an active comment that never grows stale because it evolves with the model code itself.

```julia
@mtkmodel CauerLowPassAnalog begin
    @parameters begin
        L1 = 1.304,             [description="Inductance filter coefficient no.1"]
        L2 = 0.8586,            [description="Inductance filter coefficient no.2"]
        C1 = 1.072,             [description="Capacitance filter coefficient no.1"]
        C2 = 1/(1.704992^2*L1), [description="Capacitance filter coefficient no.2"]
        C3 = 1.682,             [description="Capacitance filter coefficient no.3"]
        C4 = 1/(1.179945^2*L2), [description="Capacitance filter coefficient no.4"]
        C5 = 0.7262,            [description="Capacitance filter coefficient no.5"]
    end
    @components begin #= ... =# end
    @equations begin #= ... =# end
end
```

#### Defining Model Components

Defining the components is rather straight forward now using your paramters defined above.
Again, each of the components aside from [`Blocks.Step`](@ref) live in the [Electrical Components](@ref "ModelingToolkitStandardLibrary: Electrical Components") module.

```julia
@mtkmodel CauerLowPassAnalog begin
    @parameters begin #= ... =# end
    @components begin
        step = Step(height=1, start_time=1, smooth=false)
        source = Voltage()
        ground = Ground()
        resistor1 = Resistor(R=1)
        resistor2 = Resistor(R=1)
        inductor1 = Inductor(L=L1)
        inductor2 = Inductor(L=L2)
        capacitor1 = Capacitor(C=C1)
        capacitor2 = Capacitor(C=C2)
        capacitor3 = Capacitor(C=C3)
        capacitor4 = Capacitor(C=C4)
        capacitor5 = Capacitor(C=C5)
    end
    @equations begin #= ... =# end
end
```

#### Defining Model Equations

Defining the equations simply requires connecting all the components.
This is done with the `connect` method and knowledge of the component ports.
If you are unsure what ports are available, see the components help section named "Connectors" to learn more.

```julia
@mtkmodel CauerLowPassAnalog begin
    @parameters begin #= ... =# end
    @components begin #= ... =# end
    @equations begin
        connect(step.output, source.V)
        connect(resistor1.n, capacitor1.p)
        connect(capacitor1.n, ground.g)
        connect(inductor1.p, capacitor2.p)
        connect(inductor1.p, capacitor1.p)
        connect(inductor1.n, capacitor2.n)
        connect(capacitor2.n, capacitor3.p)
        connect(capacitor2.n, capacitor4.p)
        connect(capacitor2.n, inductor2.p)
        connect(inductor2.n, capacitor4.n)
        connect(capacitor4.n, capacitor5.p)
        connect(capacitor4.n, resistor2.p)
        connect(capacitor1.n, capacitor3.n)
        connect(capacitor1.n, capacitor5.n)
        connect(resistor2.n, capacitor1.n)
        connect(resistor1.p, source.p)
        connect(source.n, ground.g)
    end
end
```

### Defining the Problem

Now the convenience of `@mtkmodel` shines.
There is no need to declare an `ODESystem` with a `name`, a list of connections, variables and systems.
The `@mtkmodel` macro has captured all that information for you, and its complement, `@mtkbuild`, will handle the rest.

```julia
@mtkbuild model = CauerLowPassAnalog()
```

You now have a `model` which can be used to construct an `ODEProblem`.

```julia
prob = ODEProblem(model, ModelingToolkit.missing_variable_defaults(model), (0.0, 60.0))
```

What is happening with `ModelingToolkit.missing_variable_defaults(model)`?
If you try and run `ODEProblem(model, (0.0, 60.0))` then you will encounter the error:

```plaintext
ERROR: ArgumentError: Equations (7), states (7), and initial conditions (2) are of different lengths.
```

This occurs because dummy derivatives where generated without defined initial values.
Thankfully, `ModelingToolkit` provides `missing_variable_defaults` as a solution to this problem.
See [ModelingToolkit FAQ](https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/) for more information on this and other common questions.

### Solving the Problem

You are finally ready to solve!
However, if you try to run `solve(prob)` then you will encounter an error.

```plaintext
ERROR: Default algorithm choices require DifferentialEquations.jl.
Please specify an algorithm (e.g., `solve(prob, Tsit5())` or
`init(prob, Tsit5())` for an ODE) or import DifferentialEquations
directly.

You can find the list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
and its associated pages.
```

Thankfully, the error message provides instructions for next steps along with a helpful link.
The issue is resolved by explicitly specifying an algorithm, here `Rosenbrock23`.

You should now have a solution object with a successful return code.

```@example cauer_low_pass_analog
SciMLBase.successful_retcode(sol)
```
