module FluxTubes
using ModelingToolkit

@parameters t
D = Differential(t)

include("basic.jl")
include("sources.jl")

export MagneticPort, 
    Ground,
    Idle,
    Short,
    Crossing,
    ConstantPermeance,
    ConstantReluctance,
    ConstantMagneticPotentialDifference,
    ConstantMagneticFlux,
    # FluxTubes sensors
    MagneticFluxSensor, MagneticPotentialDifferenceSensor
    
end  #module