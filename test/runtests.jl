using SafeTestsets

@time @safetestset "Quality Assurance" include("qa.jl")

# Blocks
@time @safetestset "Blocks: math" include("Blocks/math.jl")
@time @safetestset "Blocks: nonlinear" include("Blocks/nonlinear.jl")
@time @safetestset "Blocks: continuous" include("Blocks/continuous.jl")
@time @safetestset "Blocks: sources" include("Blocks/sources.jl")
@time @safetestset "Blocks: analysis points" include("Blocks/test_analysis_points.jl")

# Electrical
@time @safetestset "Analog Circuits" include("Electrical/analog.jl")
@time @safetestset "Chua Circuit Demo" include("chua_circuit.jl")
@time @safetestset "Digital Circuits" include("Electrical/digital.jl")
@time @safetestset "Chua Circuit Demo" include("chua_circuit.jl")

# Thermal
@time @safetestset "Thermal Circuits" include("Thermal/thermal.jl")
@time @safetestset "Thermal Demo" include("Thermal/demo.jl")

# Magnetic
@time @safetestset "Magnetic" include("Magnetic/magnetic.jl")

# Mechanical
@time @safetestset "Mechanical Rotation" include("Mechanical/rotational.jl")
@time @safetestset "Mechanical Translation" include("Mechanical/translational.jl")
@time @safetestset "Mechanical Translation" include("Mechanical/translational_modelica.jl")
@time @safetestset "Multi-Domain" include("multi_domain.jl")

# Hydraulic
@time @safetestset "Hydraulic IsothermalCompressible" include("Hydraulic/isothermal_compressible.jl")
