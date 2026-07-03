module FluxTubes
using ModelingToolkitBase: ModelingToolkitBase, @component, @connector, @named,
    @parameters, @unpack, @variables, Equation, Flow, System, connect, extend,
    t_nounits as t, D_nounits as D
using ...Electrical: Pin

export PositiveMagneticPort, NegativeMagneticPort, TwoPort
include("utils.jl")

export Ground, Idle, Short, Crossing, ConstantPermeance, ConstantReluctance, EddyCurrent,
    ElectroMagneticConverter
include("basic.jl")

export ConstantMagneticPotentialDifference, ConstantMagneticFlux
include("sources.jl")

end  #module
