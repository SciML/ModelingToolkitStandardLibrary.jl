using Documenter, ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Mechanical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Magnetic
using ModelingToolkitStandardLibrary.Magnetic.FluxTubes
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal

makedocs(
    sitename="ModelingToolkitStandardLibrary.jl",
    authors="Julia Computing",
    clean=true,
    doctest=false,
    modules=[ModelingToolkitStandardLibrary,
             ModelingToolkitStandardLibrary.Blocks,
             ModelingToolkitStandardLibrary.Mechanical,
             ModelingToolkitStandardLibrary.Mechanical.Rotational,
             ModelingToolkitStandardLibrary.Magnetic,
             ModelingToolkitStandardLibrary.Magnetic.FluxTubes,
             ModelingToolkitStandardLibrary.Electrical,
             ModelingToolkitStandardLibrary.Thermal],

    format=Documenter.HTML(assets=["assets/favicon.ico"],
                           canonical="https://mtkstdlib.sciml.ai/stable/"),

    pages=[
        "ModelingToolkitStandardLibrary.jl: A Standard Library for ModelingToolkit" => "index.md",

        "Tutorials" => [
            "RC Circuit" => "tutorials/rc_circuit.md"
        ],

        "API" => [
            "Basic Blocks" => "API/blocks.md",
            "Electrical Components" => "API/electrical.md",
            "Magnetic Components" => "API/magnetic.md",
            "Mechanical Components" => "API/mechanical.md",
            "Thermal Components" => "API/thermal.md"
        ],
    ]
)

deploydocs(
    repo="github.com/SciML/ModelingToolkitStandardLibrary.jl";
    push_preview=true
)
