module Magnetic

using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq

@parameters t
D = Differential(t)

include("utils.jl")

include("FluxTubes/sensors.jl")

include("QuasiStatic/FluxTubes/basic.jl")
include("QuasiStatic/FluxTubes/sources.jl")
include("QuasiStatic/FluxTubes/sensors.jl")

export MagneticPort, 
    PositiveMagneticPort, 
    NegativeMagneticPort,
    Ground,
    TwoPortElementary,
    TwoPortExtended,
    TwoPort,
    Idle,
    Short,
    Crossing,
    ConstantPermeance,
    ConstantReluctance,
    ConstantMagneticPotentialDifference,
    ConstantMagneticFlux,
    AbsoluteSensor,
    # FluxTubes sensors
    MagneticFluxSensor, MagneticPotentialDifferenceSensor

end #module