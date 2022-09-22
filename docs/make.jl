using Documenter

include("pages.jl")

makedocs(sitename = "ModelingToolkitStandardLibrary.jl",
         authors = "Julia Computing",
         clean = true,
         doctest = false,
         strict = false,
         format = Documenter.HTML(assets = ["assets/favicon.ico"],
                                  canonical = "https://mtkstdlib.sciml.ai/stable/",
                                  prettyurls = false),
         pages = pages)

makedocs(sitename = "ModelingToolkitStandardLibrary.jl",
         authors = "Julia Computing",
         clean = true,
         doctest = false,
         strict = false,
         modules = [ModelingToolkitStandardLibrary,
             ModelingToolkitStandardLibrary.Blocks,
             ModelingToolkitStandardLibrary.Mechanical,
             ModelingToolkitStandardLibrary.Mechanical.Rotational,
             ModelingToolkitStandardLibrary.Magnetic,
             ModelingToolkitStandardLibrary.Magnetic.FluxTubes,
             ModelingToolkitStandardLibrary.Electrical,
             ModelingToolkitStandardLibrary.Thermal],
         format = Documenter.LaTeX(),
         pages = pages)

deploydocs(repo = "github.com/SciML/ModelingToolkitStandardLibrary.jl";
           push_preview = true)
