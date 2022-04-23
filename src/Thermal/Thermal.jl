"""
Library of thermal system components to model heat transfer.
"""
module Thermal
using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq

@parameters t
D = Differential(t)

include("utils.jl")

# Library of 1-dimensional heat transfer with lumped elements
include("HeatTransfer/ideal_components.jl")
include("HeatTransfer/sensors.jl")
include("HeatTransfer/sources.jl")

# Simple components for 1-dimensional incompressible thermo-fluid flow models
# TODO:
# - FluidHeatFlow

export # Interface
       HeatPort,
       # Thermal Components
       BodyRadiation, ConvectiveConductor, ConvectiveResistor,
       HeatCapacitor, ThermalConductor, ThermalResistor, ThermalCollector,
       # Thermal Sensors
       RelativeTemperatureSensor, HeatFlowSensor, TemperatureSensor,
       # Thermal Sources
       FixedHeatFlow, FixedTemperature, ThermalGround 

end