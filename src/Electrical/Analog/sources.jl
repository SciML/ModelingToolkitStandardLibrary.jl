"""
    Voltage(;name)

Acts as an ideal voltage source with no internal resistance.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
  - `V` [RealInput](@ref) Input for the voltage control signal, i.e. `V ~ p.v - n.v`
"""
@component function Voltage(; name)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @named V = RealInput()
    eqs = [
        v ~ V.u,
    ]

    extend(ODESystem(eqs, t, [], []; name = name, systems = [V]), oneport)
end

"""
    Current(;name)

Acts as an ideal current source with no internal resistance.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
  - `I` [RealInput](@ref) Input for the current control signal, i.e. `I ~ p.i
"""
@component function Current(; name)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @named I = RealInput()
    eqs = [
        i ~ I.u,
    ]

    extend(ODESystem(eqs, t, [], []; name = name, systems = [I]), oneport)
end
