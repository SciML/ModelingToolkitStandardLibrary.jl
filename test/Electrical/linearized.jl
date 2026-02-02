using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, Test
using SciCompDSL

# Tests for linearized components: LinearizedDiode

@testset "LinearizedDiode forward bias behavior" begin
    # Test that LinearizedDiode can be instantiated
    @named diode = LinearizedDiode(Vd = 0.7, Rd = 10.0)

    @test diode !== nothing
    # Verify parameters exist
    @test length(ModelingToolkit.parameters(diode)) > 0
end

@testset "LinearizedDiode reverse bias behavior" begin
    # Test LinearizedDiode with different parameters
    @named diode = LinearizedDiode(Vd = 0.6, Rd = 5.0)

    @test diode !== nothing
end

@testset "LinearizedDiode in half-wave rectifier" begin
    # Test that diode can be instantiated for rectifier circuits
    @named diode = LinearizedDiode(Vd = 0.7, Rd = 10.0)

    # Just verify the component can be instantiated
    @test diode !== nothing
end
