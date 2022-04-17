module FluxTubes
using ModelingToolkit
using ..Electrical: Pin

@parameters t
D = Differential(t)

export PositiveMagneticPort, NegativeMagneticPort, TwoPort
include("utils.jl")

export Ground, Idle, Short, Crossing, ConstantPermeance, ConstantReluctance, EddyCurrent, ElectroMagneticConverter
include("basic.jl")

export  ConstantMagneticPotentialDifference, ConstantMagneticFlux
include("sources.jl")

end  #module