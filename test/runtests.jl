using SafeTestsets

@safetestset "Blocks: math" begin include("test_math.jl") end
@safetestset "Blocks: function generators" begin include("wave_functions.jl") end
@safetestset "Blocks: nonlinear" begin include("test_nonlinear.jl") end
@safetestset "Blocks: continuous" begin include("test_continuous.jl") end
@safetestset "Analog Circuits" begin include("analog.jl") end
# @safetestset "Digital Circuits" begin include("digital.jl") end
@safetestset "RC Circuit Demo" begin include("demo.jl") end
@safetestset "Thermal Circuits" begin include("thermal.jl") end
