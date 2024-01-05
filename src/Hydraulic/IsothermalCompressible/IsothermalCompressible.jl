"""
Library to model iso-thermal compressible liquid fluid flow 
"""
module IsothermalCompressible

using ModelingToolkit, Symbolics

using ...Blocks: RealInput, RealOutput
using ...Mechanical.Translational: MechanicalPort, Mass

using IfElse: ifelse

@parameters t
D = Differential(t)

export HydraulicPort, HydraulicFluid
include("utils.jl")

export Cap, Tube, FixedVolume, Volume, DynamicVolume
include("components.jl")

export MassFlow, Pressure, FixedPressure
include("sources.jl")

end
