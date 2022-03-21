module Thermal
using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq

@parameters t
D = Differential(t)

include("HeatTransfer/ideal_components.jl")
include("HeatTransfer/sensors.jl")
include("HeatTransfer/sources.jl")
include("utils.jl")

    
export # Thermal Components
       BodyRadiation, ConvectiveConductor, ConvectiveResistor,
       HeatCapacitor, ThermalConductor, ThermalResistor, ThermalCollector,
       # Thermal Sensors
       RelativeTemperatureSensor, HeatFlowSensor, TemperatureSensor,
       # Thermal Sources
       FixedHeatFlow, FixedTemperature, ThermalGround,
       # Interface
       HeatPort

end