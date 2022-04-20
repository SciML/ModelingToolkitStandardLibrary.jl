"""
Generate constant signal.

# Parameters:
- `k`: Constant output value
"""
function Constant(;name, k=1)
    @named output = RealOutput()
    pars = @parameters k=k
    eqs = [
        output.u ~ k
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate sine signal.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: Amplitude of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < star_time`
"""
function Sine(;name, 
    frequency=1, 
    amplitude=1,
    phase=0,
    offset=0,
    start_time=0)

    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase
    eqs = [
        output.u ~ offset + ifelse(t < start_time, 0, amplitude* sin(2*pi*frequency*(t - start_time) + phase))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate cosine signal.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: Amplitude of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < star_time`
"""
function Cosine(;name, 
    frequency=1, 
    amplitude=1,
    phase=0,
    offset=0,
    start_time=0)

    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase
    eqs = [
        output.u ~ offset + ifelse(t < start_time, 0, amplitude* cos(2*pi*frequency*(t - start_time) + phase))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate clock signal.

# Parameters: 
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
"""
function Clock(;name, 
    offset=0, # Offset of output signal
    start_time=0)

    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time
    eqs = [
        output.u ~ offset + ifelse(t < start_time, 0, t - start_time)
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate ramp signal.

# Parameters:
- `height`: Height of ramp
- `duration`: [s] Duration of ramp (= 0.0 gives a Step)
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < star_time`
"""
function Ramp(;name, 
    offset=0,
    height=1,
    duration=1, 
    start_time=0)

    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time height=height duration=duration
    eqs = [
        output.u ~ offset + ifelse(t < start_time, 0, 
            ifelse(t < (start_time + duration), (t - start_time) * height / duration, height))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate step signal.

# Parameters:
- `height`: Height of step
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < star_time`
"""
function Step(;name, 
    offset=0, # Offset of output signal
    height=1,
    start_time=0)

    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time height=height
    eqs = [
        output.u ~ offset + ifelse(t < start_time, 0, height)
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate exponentially damped sine signal.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: Amplitude of sine wave
- `damping`: [1/s] Damping coefficient of sine wave
- `phase`: [rad] Phase of sine wave 
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < star_time`
"""
function ExpSine(;name, 
    frequency=1, 
    amplitude=1,
    damping=0.1,
    phase=0,
    offset=0,
    start_time=0)

    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase damping=damping
    eqs = [
        output.u ~ offset + ifelse(t < start_time, 0, amplitude * exp(-damping * (t - start_time)) * sin(2*pi*frequency*(t - start_time) + phase))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

# TODO:
# - Exponentials    Generate a rising and falling exponential signal
# - Pulse   Generate pulse signal of type Real
# - SawTooth    Generate saw tooth signal
# - Trapezoid   Generate trapezoidal signal of type Real