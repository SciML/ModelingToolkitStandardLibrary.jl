using SafeTestsets

# Blocks
@safetestset "Blocks: math" begin include("Blocks/math.jl") end
@safetestset "Blocks: nonlinear" begin include("Blocks/nonlinear.jl") end
@safetestset "Blocks: continuous" begin include("Blocks/continuous.jl") end
@safetestset "Blocks: sources" begin include("Blocks/sources.jl") end

# Electrical
@safetestset "Analog Circuits" begin include("Electrical/analog.jl") end
#@safetestset "Digital Circuits" begin include("Electrical/digital.jl") end
@safetestset "RC Circuit Demo" begin include("demo.jl") end

# Thermal
@safetestset "Thermal Circuits" begin include("Thermal/thermal.jl") end
@safetestset "Thermal Demo" begin include("Thermal/demo.jl") end

# Magnetic
@safetestset "Magnetic" begin include("Magnetic/magnetic.jl") end

# Mechanical
@safetestset "Mechanical" begin include("Mechanical/rotational.jl") end

# Hydraulic
@safetestset "Hydraulic" begin include("Hydraulic/hydraulic.jl") end
