using SafeTestsets

@safetestset "Blocks: math" begin include("Blocks/math.jl") end
@safetestset "Blocks: nonlinear" begin include("Blocks/nonlinear.jl") end
@safetestset "Blocks: continuous" begin include("Blocks/continuous.jl") end
@safetestset "Analog Circuits" begin include("Electrical/analog.jl") end
#@safetestset "Digital Circuits" begin include("Electrical/digital.jl") end
@safetestset "RC Circuit Demo" begin include("demo.jl") end
@safetestset "Thermal Circuits" begin include("Thermal/thermal.jl") end
