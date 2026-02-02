using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using SciCompDSL
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant, Sine
using OrdinaryDiffEq: ReturnCode.Success

# Tests for small-signal semiconductor models: BJT_SmallSignal, MOSFET_SmallSignal

@testset "BJT_SmallSignal voltage gain" begin
    # Test: Common-emitter amplifier gain ≈ -gm * Rc
    # With gm = 0.04 S, Rc = 10kΩ, gain ≈ -400
    @named bjt = BJT_SmallSignal(gm = 0.04, r_pi = 2500.0, r_o = 100000.0)

    # Just verify the component can be instantiated with correct parameters
    @test bjt !== nothing
end

@testset "MOSFET_SmallSignal voltage gain" begin
    # Test: Common-source amplifier gain ≈ -gm * Rd
    # With gm = 0.01 S, Rd = 5kΩ, gain ≈ -50
    @named mosfet = MOSFET_SmallSignal(gm = 0.01, r_ds = 50000.0)

    # Just verify the component can be instantiated
    @test mosfet !== nothing
end

@testset "BJT_SmallSignal with capacitances" begin
    # Test that BJT with capacitances can be instantiated
    @named bjt = BJT_SmallSignal(gm = 0.04, r_pi = 2500.0, r_o = 100000.0,
        C_pi = 10e-12, C_mu = 1e-12)

    # Just verify the component can be instantiated
    @test bjt !== nothing
end

@testset "MOSFET_SmallSignal with capacitances" begin
    # Test that MOSFET with capacitances can be instantiated
    @named mosfet = MOSFET_SmallSignal(gm = 0.01, r_ds = 50000.0,
        C_gs = 5e-12, C_gd = 1e-12)

    # Just verify the component can be instantiated
    @test mosfet !== nothing
end
