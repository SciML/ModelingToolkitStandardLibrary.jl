using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit

@testset "Guesses" begin
    for (block, guess) in [
            (RealInput(; name = :a), 0.0),
            (RealInput(; nin = 3, name = :a), zeros(3)),
            (RealOutput(; name = :a), 0.0),
            (RealOutput(; nout = 3, name = :a), zeros(3)),
            (RealInputArray(; nin = 3, name = :a), zeros(3)),
            (RealOutputArray(; nout = 3, name = :a), zeros(3)),
        ]
        guesses = ModelingToolkit.guesses(block)
        @test guesses[@nonamespace block.u] == guess
    end
end
