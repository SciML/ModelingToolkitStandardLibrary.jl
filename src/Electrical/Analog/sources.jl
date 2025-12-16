"""
    Voltage(; name)

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

    pars = @parameters begin
    end

    systems = @named begin
        V = RealInput()
    end

    vars = @variables begin
    end

    equations = Equation[
        v ~ V.u
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
end

"""
    Current(; name)

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

    pars = @parameters begin
    end

    systems = @named begin
        I = RealInput()
    end

    vars = @variables begin
    end

    equations = Equation[
        i ~ I.u
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
end
