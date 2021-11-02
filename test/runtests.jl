using SafeTestsets

@safetestset "Analog Circuits" begin include("analog.jl") end
#@safetestset "Digital Circuits" begin include("digital.jl") end
@safetestset "RC Circuit Demo" begin include("demo.jl") end
@safetestset "Thermal Circuits" begin include("thermal.jl") end
