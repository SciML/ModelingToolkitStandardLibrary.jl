module FluxTubes
using ModelingToolkit
using ModelingToolkit: t, D
using ...Electrical: Pin
import ...Wb
using ...DynamicQuantities: @u_str

export PositiveMagneticPort, NegativeMagneticPort, TwoPort
include("utils.jl")

export Ground, Idle, Short, Crossing, ConstantPermeance, ConstantReluctance, EddyCurrent,
       ElectroMagneticConverter
include("basic.jl")

export ConstantMagneticPotentialDifference, ConstantMagneticFlux
include("sources.jl")

end  #module
