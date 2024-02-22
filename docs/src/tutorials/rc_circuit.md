# RC Circuit Model

This tutorial is a simplified version of the [RC circuit tutorial in the
`ModelingToolkit.jl` documentation](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/acausal_components/).
In that tutorial, the full RC circuit is built from scratch. Here, we will use the
components of the `Electrical` model in the ModelingToolkit Standard Library to simply
connect pre-made components and simulate the model.

```@example
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant

R = 1.0
C = 1.0
V = 1.0
@variables t
@named resistor = Resistor(R = R)
@named capacitor = Capacitor(C = C, v = 0.0)
@named source = Voltage()
@named constant = Constant(k = V)
@named ground = Ground()

rc_eqs = [connect(constant.output, source.V)
          connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n, ground.g)]

@named rc_model = ODESystem(rc_eqs, t,
    systems = [resistor, capacitor, constant, source, ground])
sys = structural_simplify(rc_model)
prob = ODEProblem(sys, Pair[], (0, 10.0))
sol = solve(prob, Tsit5())
plot(sol, idxs = [capacitor.v, resistor.i],
    title = "RC Circuit Demonstration",
    labels = ["Capacitor Voltage" "Resistor Current"])
savefig("plot.png");
nothing; # hide
```

![](plot.png)
