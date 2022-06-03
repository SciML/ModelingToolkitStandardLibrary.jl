"""
Library of electrical models.
This library contains electrical components to build up analog circuits.
"""
module Electrical

using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq
using OffsetArrays

@parameters t
D = Differential(t)

using ..Blocks: RealInput, RealOutput

include("utils.jl")
include("Analog/ideal_components.jl")
include("Analog/sensors.jl")
include("Analog/sources.jl")
# include("Digital/components.jl")
# include("Digital/gates.jl")
# include("Digital/tables.jl")
# include("Digital/sources.jl")

# TODO:
# - digital
# - machines
# - multi-phase

export #Interface
       Pin,
       # Analog Components
       Capacitor, Ground, Inductor, Resistor, Conductor,
       Short, IdealOpAmp,
       # Analog Sensors
       CurrentSensor, PotentialSensor, VoltageSensor,
       PowerSensor, MultiSensor,
       # Analog Sources
       Voltage, Current


       # # Digital Gates
       # And, Or, Not, Xor, Nand, Nor, Xnor,
       # # Digital components
       # HalfAdder, FullAdder, MUX, DEMUX, Encoder, Decoder,
       # # Digital Sources
       # DigitalPin, Pulse, PulseDiff

end
