using ModelingToolkitStandardLibrary, JET, Test

# JET.jl static analysis tests for ModelingToolkitStandardLibrary
#
# Note: This package heavily uses ModelingToolkit's @component and @connector macros
# which generate code at macro expansion time. Full package analysis with JET times out.
# Instead, we focus on testing specific utility functions that are not macro-generated.

@testset "JET static analysis" begin
    # Test that the package loads without JET-detectable issues
    # We use a targeted approach since full package analysis is not feasible
    # for metaprogramming-heavy packages like this one.

    @testset "Module loading" begin
        # Verify modules can be accessed (basic sanity check)
        @test isdefined(ModelingToolkitStandardLibrary, :Blocks)
        @test isdefined(ModelingToolkitStandardLibrary, :Electrical)
        @test isdefined(ModelingToolkitStandardLibrary, :Mechanical)
        @test isdefined(ModelingToolkitStandardLibrary, :Thermal)
        @test isdefined(ModelingToolkitStandardLibrary, :Magnetic)
        @test isdefined(ModelingToolkitStandardLibrary, :Hydraulic)
    end

    @testset "Utility functions - Hydraulic" begin
        # Test hydraulic utility functions for type stability
        # These are pure numerical functions that don't use metaprogramming

        # regPow function - regularized power
        regPow = ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible.regPow
        @test_opt target_modules = (ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible,) regPow(
            1.0, 0.5, 0.01
        )

        # transition function - smooth transition between values
        transition = ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible.transition
        @test_opt target_modules = (ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible,) transition(
            0.0, 1.0, 0.0, 1.0, 0.5
        )

        # friction_factor function - calculates friction factor for pipe flow
        friction_factor = ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible.friction_factor
        @test_opt target_modules = (ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible,) friction_factor(
            1.0, 0.01, 0.1, 0.001, 64.0
        )
    end

    @testset "Utility functions - Blocks sources" begin
        # Test smooth functions used in source blocks
        # These are pure numerical functions for generating smooth waveforms

        Blocks = ModelingToolkitStandardLibrary.Blocks

        # smooth_step - smooth step function
        @test_opt target_modules = (Blocks,) Blocks.smooth_step(1.0, 0.01, 1.0, 0.0, 0.0)

        # smooth_xH - smooth ramp helper
        @test_opt target_modules = (Blocks,) Blocks.smooth_xH(1.0, 0.01, 0.0)

        # square - square wave
        @test_opt target_modules = (Blocks,) Blocks.square(1.0, 1.0, 1.0, 0.0, 0.0)

        # triangular - triangular wave
        @test_opt target_modules = (Blocks,) Blocks.triangular(1.0, 1.0, 1.0, 0.0, 0.0)
    end

    @testset "Digital logic - Electrical" begin
        # Test digital logic table operations
        E = ModelingToolkitStandardLibrary.Electrical

        # Logic enum operations
        @test E.F0 isa E.Logic
        @test E.F1 isa E.Logic

        # Test logic operations type stability
        @test_opt target_modules = (E,) E._and2(E.F0, E.F1)
        @test_opt target_modules = (E,) E._or2(E.F0, E.F1)
        @test_opt target_modules = (E,) E._not(E.F0)
        @test_opt target_modules = (E,) E._xor2(E.F0, E.F1)
    end
end
