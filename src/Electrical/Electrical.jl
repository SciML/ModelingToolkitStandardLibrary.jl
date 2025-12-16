"""
Library of electrical models.
This library contains electrical components to build up analog circuits.
"""
module Electrical

using ModelingToolkit, Symbolics, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using ..Thermal: HeatPort
using ..Mechanical.Rotational: Flange, Support
using ..Blocks: RealInput, RealOutput

export Pin, OnePort, DigitalPin
include("utils.jl")

export Capacitor,
       Ground, Inductor, Resistor, Conductor, Short, IdealOpAmp, EMF,
       Diode, VariableResistor
include("Analog/ideal_components.jl")

export CurrentSensor, PotentialSensor, VoltageSensor, PowerSensor, MultiSensor
include("Analog/sensors.jl")

export Voltage, Current
include("Analog/sources.jl")

export NMOS, PMOS
include("Analog/mosfets.jl")

export NPN, PNP
include("Analog/transistors.jl")

export Logic
include("Digital/logic.jl")

export StdLogicVector, StdULogicVector,
       std_ulogic, UX01, UX01Z, X01, X01Z,
       get_logic_level
include("Digital/logic_vectors.jl")

export LogicTable,
       AndTable, OrTable, NotTable, XorTable,
       _not, _and, _or, _xor
include("Digital/tables.jl")

export Not, And, Nand, Or, Nor, Xor, Xnor
include("Digital/gates.jl")

export HalfAdder, FullAdder, MUX, DEMUX, Encoder, Decoder
include("Digital/components.jl")

export PulseDiff, Set, Reset, Pulse
include("Digital/sources.jl")

# TODO:
# - digital
# - machines
# - multi-phase

end
