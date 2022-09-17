"""
Library of mechanical models.
"""
module Mechanical

using ModelingToolkit

include("Rotational/Rotational.jl")
include("Translational/Translational.jl")
include("TranslationalPosition/TranslationalPosition.jl")

end
