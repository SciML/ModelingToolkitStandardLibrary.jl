# ModelingToolkitStandardLibrary.jl

[![CI](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/actions/workflows/CI.yml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://mtkstdlib.sciml.ai/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://mtkstdlib.sciml.ai/dev/)

The ModelingToolkit Standard Library is a standard library of components to model the world and beyond.

![](https://user-images.githubusercontent.com/1814174/172000112-3579f5cf-c370-48c2-8047-558fbc46aeb6.png)

## Installation

Assuming that you already have Julia correctly installed, it suffices to import
ModelingToolkitStandardLibrary.jl in the standard way:

```julia
import Pkg; Pkg.add("ModelingToolkitStandardLibrary")
```

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://mtkstdlib.sciml.ai/stable/). Use the
[in-development documentation](https://mtkstdlib.sciml.ai/dev/) for the version of
the documentation, which contains the unreleased features.

## Libraries

The following are the constituant libraries of the ModelingToolkit Standard Library.

- [Basic Blocks](http://mtkstdlib.sciml.ai/dev/API/blocks/)
- [Mechanical Components](http://mtkstdlib.sciml.ai/dev/API/mechanical/)
- [Electrical Components](http://mtkstdlib.sciml.ai/dev/API/electrical/)
- [Magnetic Components](http://mtkstdlib.sciml.ai/dev/API/magnetic/)
- [Thermal Components](http://mtkstdlib.sciml.ai/dev/API/thermal/)

## Example

The following is the [RC Circuit Demonstration](http://mtkstdlib.sciml.ai/dev/tutorials/rc_circuit/):

```julia
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant

R = 1.0
C = 1.0
V = 1.0
@variables t
@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = Voltage()
@named constant = Constant(k=V)
@named ground = Ground()

rc_eqs = [
        connect(constant.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
        ]

@named rc_model = ODESystem(rc_eqs, t, systems=[resistor, capacitor, constant, source, ground])
sys = structural_simplify(rc_model)
prob = ODAEProblem(sys, Pair[], (0, 10.0))
sol = solve(prob, Tsit5())
plot(sol, vars = [capacitor.v, resistor.i],
     title = "RC Circuit Demonstration",
     labels = ["Capacitor Voltage" "Resistor Current"])
savefig("plot.png")
```

![](https://user-images.githubusercontent.com/1814174/164912983-c3f73628-0e19-4e42-b085-4f62ba6f23d1.png)
