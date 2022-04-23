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
    ConstantVoltage(;name, V = 1.0)   

Source for constant voltage,

# Parameters:
- `V`: [V] Voltage
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
    CosineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

Generate cosine voltage.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: [V] Amplitude of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: [V] Offset of output voltage
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
Generate damped sine voltage.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: [V] Amplitude of sine wave
- `damping`: [1/s] Damping coefficient of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: [V] Offset of output voltage
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    RampVoltage(;name, offset=0.0, start_time=0.0, duration=1.0, height=1.0)

Generate ramp voltage.

# Parameters:
- `height`: [V] Height of ramp
- `duration`: [s] Duration of ramp (= 0.0 gives a Step)
- `offset`: [V] Offset of output voltage
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    SineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

Generate sine voltage.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: [V] Amplitude of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: [V] Offset of output voltage
- `start_time`: [s] Output `y = offset` for `t < start_time`
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

Generate square voltage.
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
    StepVoltage(;name, offset=0.0, start_time=0.0, height=1.0)

Generate step voltage.

# Parameters:
- `height`: [V] Height of step
- `offset`: [V] Offset of output voltage
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    ConstantCurrent(;name, I = 1.0)   

Source for constant current.

# Parameters:
- `I`: [A] Current
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
    CosineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

Generate cosine current.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: [A] Amplitude of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: [A] Offset of output current
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    ExpSineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0, damping=0.0)

Generate damped sine current.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: [A] Amplitude of sine wave
- `damping`: [1/s] Damping coefficient of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: [A] Offset of output current
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    RampCurrent(;name, offset=0.0, start_time=0.0, duration=1.0, height=1.0)

Generate ramp current.

# Parameters:
- `height`: [A] Height of ramp
- `duration`: [s] Duration of ramp (= 0.0 gives a Step)
- `offset`: [A] Offset of output current
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    SineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, start_time=0.0, phase=0.0)

Generate sine current.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: [A] Amplitude of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: [A] Offset of output current
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
    StepCurrent(;name, offset=0.0, start_time=0.0, height=1.0)

Generate step current.

# Parameters:
- `height`: [A] Height of step
- `offset`: [A] Offset of output current
- `start_time`: [s] Output `y = offset` for `t < start_time`
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
