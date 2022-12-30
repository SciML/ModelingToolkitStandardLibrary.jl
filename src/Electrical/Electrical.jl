"""
Library of electrical models.
This library contains electrical components to build up analog circuits.
"""
module Electrical

using ModelingToolkit, Symbolics, IfElse
using ..Thermal: HeatPort
using ..Mechanical.Rotational: Flange, Support
using ..Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Pin, OnePort
include("utils.jl")

export Capacitor, Ground, Inductor, Resistor, Conductor, Short, IdealOpAmp, EMF,
       HeatingResistor
include("Analog/ideal_components.jl")

export CurrentSensor, PotentialSensor, VoltageSensor, PowerSensor, MultiSensor
include("Analog/sensors.jl")

export Voltage, Current
include("Analog/sources.jl")

# include("Digital/components.jl")
# include("Digital/gates.jl")
# include("Digital/sources.jl")

# TODO:
# - digital
# - machines
# - multi-phase

export Logic
include("Digital/logic.jl")

include("Digital/tables.jl")

end
