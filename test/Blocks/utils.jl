using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit

@testset "Guesses" begin
    for (block, guess) in [
        (RealInputArray(; nin = 3, name = :a), zeros(3)),
        (RealOutputArray(; nout = 3, name = :a), zeros(3))
    ]
        guesses = ModelingToolkit.guesses(block)
        @test guesses[@nonamespace block.u] == guess
    end
end

@testset "Guesses" begin
    for (block, guess) in [
        (RealInput(; name = :a), 0.0),
        (RealInput(; nin = 3, name = :a), zeros(3)),
        (RealOutput(; name = :a), 0.0),
        (RealOutput(; nout = 3, name = :a), zeros(3)),
    ]
        guesses = ModelingToolkit.guesses(block)
        @test guesses[@nonamespace block.u[1]] == guess[1]
    end
end

@test_deprecated RealInput(; name = :a, u_start = 1.0)
@test_deprecated RealInput(; name = :a, nin = 2, u_start = ones(2))
@test_deprecated RealOutput(; name = :a, u_start = 1.0)
@test_deprecated RealOutput(; name = :a, nout = 2, u_start = ones(2))
@test_deprecated RealInputArray(; name = :a, nin = 2, u_start = ones(2))
@test_deprecated RealOutputArray(; name = :a, nout = 2, u_start = ones(2))
