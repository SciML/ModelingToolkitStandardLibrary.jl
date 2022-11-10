# ModelingToolkitStandardLibrary.jl

ModelingToolkitStandardLibrary.jl is a standard library for the 
[ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/) acasual modeling system.

## Installation

Assuming that you already have Julia correctly installed, it suffices to import
ModelingToolkitStandardLibrary.jl in the standard way:

```julia
import Pkg; Pkg.add("ModelingToolkitStandardLibrary")
```

## Tutorials 

- [RC Circuit](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/)
- [Custom Component](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/custom_component/)
- [Thermal Model](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/thermal_model/)
- [DC Motor with PI-controller](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/)

## Libraries

The following are the constituant libraries of the ModelingToolkit Standard Library.

- [Basic Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/)
- [Mechanical Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/mechanical/)
- [Electrical Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/electrical/)
- [Magnetic Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/magnetic/)
- [Thermal Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/thermal/)

## Reproducibility
```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```
```@example
using Pkg # hide
Pkg.status() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>and using this machine and Julia version.</summary>
```
```@example
using InteractiveUtils # hide
versioninfo() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```
```@example
using Pkg # hide
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```
```@raw html
</details>
```
```@raw html
You can also download the 
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```
```@raw html
">project</a> file.
```