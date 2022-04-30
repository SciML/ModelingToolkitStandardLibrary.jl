# Define and register smooth functions
_cos_wave(t, f, A, st, ϕ) = A*cos(2*π*f*(t - st) + ϕ)
_sin_wave(t, f, A, st, ϕ) = A*sin(2*π*f*(t - st) + ϕ)
_damped_sine_wave(t, f, A, st, ϕ, d) = exp((st-t)*d)*A*sin(2*π*f*(t-st) + ϕ)
_ramp(t, δ, st, et, h) = h/(et-st)*(_xH(t, δ, st) - _xH(t, δ, et))
_square_wave(t, δ, f, A, st) = A*2atan(sin(2π*(t-st)*f)/δ)/π
_step(t, δ, h, a) = h*(atan((t-a)/δ)/π + 0.5)
_triangular_wave(t, δ, f, A, st) = A*(1-2acos((1 - δ)sin(2π*(t-st)*f))/π)
_xH(t, δ, tₒ) = (t-tₒ)*(1+((t-tₒ)/sqrt((t-tₒ)^2+δ^2)))/2

@register_symbolic _cos_wave(t, f, A, st, ϕ)
@register_symbolic _sin_wave(t, f, A, st, ϕ)
@register_symbolic _damped_sine_wave(t, f, A, st, ϕ, damping)
@register_symbolic _ramp(t, δ, st, et, h)
@register_symbolic _square_wave(t, δ, f, A, st)
@register_symbolic _step(t, δ, h, a)
@register_symbolic _triangular_wave(t, δ, f, A, st)

"""
    ConstantVoltage(; name, V=1.0)

The source for an ideal constant voltage.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `V`: [`V`] The constant voltage across the terminals of this source
"""
function ConstantVoltage(;name, V = 1.0)   
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters V=V
    eqs = [
        v ~ V
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    CosineVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)

A source in which the voltage across its terminals is a cosine function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Observables
- `offset`: [`V`] A constant offset added to the voltage output
- `amplitude`: [`V`] The amplitude of the cosine function
- `frequency`: [`Hz`] The frequency of the cosine function
- `starttime`: [`s`] The time at which the source starts functioning. Before this time, the voltage across its terminals is 0.
- `phase`: [`rad`] The phase offset of the cosine function
"""
function CosineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        start_time=start_time 
        phase=phase
    end
    eqs = [
        v ~ _cos_wave(t, frequency, amplitude, start_time, phase) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    ExpSineVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0, damping=0.0)

A source in which the voltage across its terminals is a damped sine function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`V`] A constant offset added to the voltage output
- `amplitude`: [`V`] The amplitude of the damped sine function
- `frequency`: [`Hz`] The frequency of the damped sine function
- `start_time`: [`s`] The time at which the source starts functioning. Before this time, the voltage across its terminals is `offset`.
- `phase`: [`rad`] The phase offset of the damped sine function
- `damping_coef`: [`1/s`] Damping coefficient of the damped sine function
"""
function ExpSineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0, damping=0.0)
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency
        start_time=start_time
        phase=phase 
        damping=damping
    end
    eqs = [
        v ~ _damped_sine_wave(t, frequency, amplitude, start_time, phase, damping) * _step(t, δ, 1.0, start_time)
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    RampVoltage(; name, offset=0.0, start_time=0.0, duration=1.0, height=1.0)

A source in which the voltage across grows linearly from `offset` to `offset+height` over
the time interval `duration` starting at `start_time`

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`V`] A constant offset added to the voltage output
- `start_time`: [`s`] The time at which the voltage starts growing
- `duration`: [`s`] The duration of the ramp (`0.0` gives a step)
- `height`: [`V`] The amount that the voltage grows in the time interval
"""
function RampVoltage(;name, offset=0.0, start_time=0.0, duration=1.0, height=1.0)
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        start_time=start_time 
        duration=duration
    end
    eqs = [
        v ~ _ramp(t, δ, start_time, start_time + duration, height) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    SineVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

A source in which the voltage across its terminals is a sine function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`V`] A constant offset added to the voltage output
- `amplitude`: [`V`] The amplitude of the sine function
- `frequency`: [`Hz`] The frequency of the sine function
- `start_time`: [`s`] The time at which the source starts functioning. Before this time, the voltage across its terminals is `offset`.
- `phase`: [`rad`] The phase offset of the sine function
"""
function SineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        start_time=start_time 
        phase=phase
    end
    eqs = [
        v ~ _sin_wave(t, frequency, amplitude, start_time, phase) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    SquareVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)

A source in which the voltage across its terminals is a square function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`V`] A constant offset added to the voltage output
- `amplitude`: [`V`] The amplitude of the square wave function
- `frequency`: [`Hz`] The frequency of the square wave function
- `start_time`: [`s`] The time at which the source starts functioning. Before this time, the voltage across its terminals is `offset`.
"""
function SquareVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        start_time=start_time 
    end
    eqs = [
        v ~ _square_wave(t, δ, frequency, amplitude, start_time) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    StepVoltage(; name, offset=0.0, start_time=0.0, height=1.0)

A source in which the voltage across its terminals increases from `offset` to `offset+height` at `start_time`

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Observables
- `offset`: [`V`] A constant offset added to the voltage output
- `start_time`: [`s`] The time at which the source starts functioning, and the voltage jumps
- `height`: [`V`] Magnitude of increase in voltage
"""
function StepVoltage(;name, offset=0.0, start_time=0.0, height=1.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        start_time=start_time 
    end
    eqs = [
        v ~ _step(t, δ, height, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    TriangularVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)

A source in which the voltage across its terminals is a triangular function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Observables
- `offset`: [`V`] A constant offset added to the voltage output
- `amplitude`: [`V`] Amplitude of the triangular wave function
- `frequency`: [`Hz`] Frequency of the triangular wave function
- `start_time`: [`s`] The time at which the source starts functioning. Before this, the output of the source is `offset`
"""
function TriangularVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        start_time=start_time 
    end
    eqs = [
        v ~ _triangular_wave(t, δ, frequency, amplitude, start_time) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

# Current Sources ######################################################################################################
"""
    ConstantCurrent(; name, I = 1.0)

The source for an ideal constant current.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `I`: [`A`] The constant current through the terminals of this source
"""
function ConstantCurrent(;name, I = 1.0)   
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters I=I
    eqs = [
        i ~ I
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    CosineCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

A source in which the current through its terminals is a cosine function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Observables
- `offset`: [`A`] A constant offset added to the current output
- `amplitude`: [`A`] The amplitude of the cosine function
- `frequency`: [`Hz`] The frequency of the cosine function
- `starttime`: [`s`] The time at which the source starts functioning. Before this time, the current through its terminals is 0.
- `phase`: [`rad`] The phase offset of the cosine function
"""
function CosineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        start_time=start_time 
        phase=phase
    end
    eqs = [
        i ~ _cos_wave(t, frequency, amplitude, start_time, phase) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    ExpSineCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0, damping=0.0)

A source in which the current through its terminals is a damped sine function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`A`] A constant offset added to the current output
- `amplitude`: [`A`] The amplitude of the damped sine function
- `frequency`: [`Hz`] The frequency of the damped sine function
- `start_time`: [`s`] The time at which the source starts functioning. Before this time, the current through its terminals is `offset`.
- `phase`: [`rad`] The phase offset of the damped sine function
- `damping_coef`: [`1/s`] Damping coefficient of the damped sine function
"""
function ExpSineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0, damping=0.0)
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency
        start_time=start_time
        phase=phase 
        damping=damping
    end
    eqs = [
        i ~ _damped_sine_wave(t, frequency, amplitude, start_time, phase, damping) * _step(t, δ, 1.0, start_time)
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    RampCurrent(; name, offset=0.0, start_time=0.0, duration=1.0, height=1.0)

A source in which the current grows linearly from `offset` to `offset+height` over
the time interval `duration` starting at `start_time`

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`A`] A constant offset added to the current output
- `start_time`: [`s`] The time at which the current starts growing
- `duration`: [`s`] The duration of the ramp (`0.0` gives a step)
- `height`: [`A`] The amount that the current grows in the time interval
"""
function RampCurrent(;name, offset=0.0, start_time=0.0, duration=1.0, height=1.0)
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        start_time=start_time 
        duration=duration
    end
    eqs = [
        i ~ _ramp(t, δ, start_time, start_time + duration, height) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    SineCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

A source in which the current through its terminals is a sine function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`A`] A constant offset added to the current output
- `amplitude`: [`V`] The amplitude of the sine function
- `frequency`: [`Hz`] The frequency of the sine function
- `start_time`: [`s`] The time at which the source starts functioning. Before this time, the current through its terminals is `offset`.
- `phase`: [`rad`] The phase offset of the sine function
"""
function SineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        start_time=start_time 
        phase=phase
    end
    eqs = [
        i ~ _sin_wave(t, frequency, amplitude, start_time, phase) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    SquareCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)

A source in which the current through its terminals is a square function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Parameters
- `offset`: [`A`] A constant offset added to the current output
- `amplitude`: [`A`] The amplitude of the square wave function
- `frequency`: [`Hz`] The frequency of the square wave function
- `start_time`: [`s`] The time at which the source starts functioning. Before this time, the current through its terminals is `offset`.
"""
function SquareCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        start_time=start_time 
    end
    eqs = [
        i ~ _square_wave(t, δ, frequency, amplitude, start_time) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    StepCurrent(; name, offset=0.0, start_time=0.0, height=1.0)

A source in which the current through its terminals increases from `offset` to `offset+height` at `start_time`

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Observables
- `offset`: [`A`] A constant offset added to the current output
- `start_time`: [`s`] The time at which the source starts functioning, and the current jumps
- `height`: [`A`] Magnitude of increase in current
"""
function StepCurrent(;name, offset=0.0, start_time=0.0, height=1.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        start_time=start_time 
    end
    eqs = [
        i ~ _step(t, δ, height, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    TriangularCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)

A source in which the current through its terminals is a triangular function of time.

# States
- See [OnePort](@ref)

# Connectors
- `p` Positive port
- `n` Negative port

# Observables
- `offset`: [`A`] A constant offset added to the current output
- `amplitude`: [`A`] Amplitude of the triangular wave function
- `frequency`: [`Hz`] Frequency of the triangular wave function
- `start_time`: [`s`] The time at which the source starts functioning. Before this, the output of the source is `offset`
"""
function TriangularCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        start_time=start_time 
    end
    eqs = [
        i ~ _triangular_wave(t, δ, frequency, amplitude, start_time) * _step(t, δ, 1.0, start_time) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end
