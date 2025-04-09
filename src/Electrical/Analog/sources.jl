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
@mtkmodel Voltage begin
    @extend v, i = oneport = OnePort()
    @components begin
        V = RealInput()
    end
    @equations begin
        v ~ V.u
    end
end

"""
    ConstantVoltage(; name, V)

Acts as an ideal constant voltage source with no internal resistance.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `V`: [`V`] Constant voltage value
"""
@mtkmodel ConstantVoltage begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        V, [description = "Constant voltage value"]
    end
    @equations begin
        v ~ V
    end
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
@mtkmodel Current begin
    @extend v, i = oneport = OnePort()
    @components begin
        I = RealInput()
    end
    @equations begin
        i ~ I.u
    end
end

"""
    ConstantCurrent(; name, I)

Acts as an ideal constant current source with no internal resistance.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `I`: [`A`] Constant current value
"""
@mtkmodel ConstantCurrent begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        I, [description = "Constant current value"]
    end
    @equations begin
        i ~ I
    end
end
