# Define and register smooth functions
# These are "smooth" aka differentiable and avoid Gibbs effect
# These follow: `offset` + `smooth_wave` * `smooth_step` with zero output for `t < start_time`
function smooth_cos(x, δ, f, amplitude, ϕ, offset, start_time)
    offset +
    amplitude * cos(2 * π * f * (x - start_time) + ϕ) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_damped_sin(x, δ, f, amplitude, damping, ϕ, offset, start_time)
    offset +
    exp((start_time - x) * damping) * amplitude * sin(2 * π * f * (x - start_time) + ϕ) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_ramp(x, δ, height, duration, offset, start_time)
    offset +
    height / (duration) *
    (smooth_xH(x, δ, start_time) - smooth_xH(x, δ, start_time + duration))
end

function smooth_sin(x, δ, f, amplitude, ϕ, offset, start_time)
    offset +
    amplitude * sin(2 * pi * f * (x - start_time) + ϕ) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_square(x, δ, f, amplitude, offset, start_time)
    offset +
    amplitude * 2atan(sin(2π * (x - start_time) * f) / δ) / π *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_step(x, δ, height, offset, start_time)
    offset + height * (atan((x - start_time) / δ) / π + 0.5)
end

function smooth_triangular(x, δ, f, amplitude, offset, start_time)
    offset +
    amplitude * (1 - 2acos((1 - δ)sin(2π * (x - start_time) * f)) / π) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_xH(x, δ, tₒ)
    0.5 * (x - tₒ) * (1 + ((x - tₒ) / sqrt((x - tₒ)^2 + δ^2)))
end

function square(x, f, amplitude, offset, start_time)
    offset +
    (x > start_time) * (amplitude *
     (4 * floor(f * (x - start_time)) - 2 * floor(2 * (x - start_time) * f) + 1))
end

function triangular(x, f, amplitude, offset, start_time)
    p = 1 / f # period
    offset +
    (x > start_time) *
    (4 * amplitude * f * abs(abs((x - p / 4 - start_time) % p) - p / 2) - amplitude)
end

"""
Generate constant signal.

# Parameters:
- `k`: Constant output value

# Connectors:
- `output`
"""
function Constant(; name, k = 1)
    @named output = RealOutput()
    pars = @parameters k = k
    eqs = [
        output.u ~ k,
    ]
    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate sine signal.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: Amplitude of sine wave
- `phase`: [rad] Phase of sine wave
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""
function Sine(; name,
              frequency,
              amplitude = 1,
              phase = 0,
              offset = 0,
              start_time = 0,
              smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase
    equation = if smooth == false
        offset + ifelse(t < start_time, 0,
               amplitude * sin(2 * pi * frequency * (t - start_time) + phase))
    else
        δ = 1e-5
        smooth_sin(t, δ, frequency, amplitude, phase, offset, start_time)
    end

    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate cosine signal.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: Amplitude of sine wave
- `phase`: [rad] Phase of sine wave
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""

function Cosine(; name,
                frequency,
                amplitude = 1,
                phase = 0,
                offset = 0,
                start_time = 0,
                smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase
    equation = if smooth == false
        offset + ifelse(t < start_time, zero(t),
               amplitude * cos(2 * pi * frequency * (t - start_time) + phase))
    else
        δ = 1e-5
        smooth_cos(t, δ, frequency, amplitude, phase, offset, start_time)
    end
    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate current time signal.

# Parameters:
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`

# Connectors:
- `output`
"""
function ContinuousClock(; name, offset = 0, start_time = 0)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time
    eqs = [
        output.u ~ offset + ifelse(t < start_time, zero(t), t - start_time),
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate ramp signal.

# Parameters:
- `height`: Height of ramp
- `duration`: [s] Duration of ramp (= 0.0 gives a Step)
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""
function Ramp(; name,
              height = 1,
              duration = 1,
              offset = 0,
              start_time = 0,
              smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time height=height duration=duration
    equation = if smooth == false
        offset + ifelse(t < start_time, 0,
               ifelse(t < (start_time + duration), (t - start_time) * height / duration,
                      height))
    else
        δ = 1e-5
        smooth_ramp(t, δ, height, duration, offset, start_time)
    end

    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate smooth square signal.

# Parameters:
- `frequency`: [Hz] Frequency of square wave
- `amplitude`: Amplitude of square wave
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""
function Square(; name, frequency = 1.0, amplitude = 1.0,
                offset = 0.0, start_time = 0.0, smooth = false)
    δ = 1e-5

    @named output = RealOutput()
    pars = @parameters begin
        frequency = frequency
        amplitude = amplitude
        offset = offset
        start_time = start_time
    end

    equation = if smooth == false
        square(t, frequency, amplitude, offset, start_time)
    else
        δ = 1e-5
        smooth_square(t, δ, frequency, amplitude, offset, start_time)
    end

    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
    Step(;name, height=1, offset=0, start_time=0, duration=Inf, smooth=true)

Generate step signal.

# Parameters:
- `height`: Height of step
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time` and thereafter `offset+height`.
- `duration`: [s] If `duration < Inf` is supplied, the output will revert to `offset` after `duration` seconds.
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""
function Step(; name, height = 1, offset = 0, start_time = 0, duration = Inf, smooth = true)
    @named output = RealOutput()
    duration_numeric = duration
    pars = @parameters offset=offset start_time=start_time height=height duration=duration
    equation = if !smooth
        offset + ifelse((start_time < t) & (t < start_time + duration), height, 0)
    else
        δ = 1e-5
        if duration_numeric == Inf
            smooth_step(t, δ, height, offset, start_time)
        else
            smooth_step(t, δ, height, offset, start_time) -
            smooth_step(t, δ, height, 0, start_time + duration)
        end
    end

    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate exponentially damped sine signal.

# Parameters:
- `frequency`: [Hz] Frequency of sine wave
- `amplitude`: Amplitude of sine wave
- `damping`: [1/s] Damping coefficient of sine wave
- `phase`: [rad] Phase of sine wave
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""
function ExpSine(; name,
                 frequency,
                 amplitude = 1,
                 damping = 0.1,
                 phase = 0,
                 offset = 0,
                 start_time = 0,
                 smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase damping=damping

    equation = if smooth == false
        offset + ifelse(t < start_time, 0,
               amplitude * exp(-damping * (t - start_time)) *
               sin(2 * pi * frequency * (t - start_time) + phase))
    else
        δ = 1e-5
        smooth_damped_sin(t, δ, frequency, amplitude, damping, phase, offset, start_time)
    end

    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
Generate smooth triangular signal for frequencies less than or equal to 25 Hz

# Parameters:
- `frequency`: [Hz] Frequency of square wave
- `amplitude`: Amplitude of square wave
- `offset`: Offset of output signal.
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a smoothing factor of `δ=1e-5`

# Connectors:
- `output`
"""
function Triangular(; name, amplitude = 1.0, frequency = 1.0,
                    offset = 0.0, start_time = 0.0, smooth = false)
    @named output = RealOutput()
    pars = @parameters begin
        amplitude = amplitude
        frequency = frequency
        offset = offset
        start_time = start_time
    end

    equation = if smooth == false
        triangular(t, frequency, amplitude, offset, start_time)
    else
        δ = 1e-5
        smooth_triangular(t, δ, frequency, amplitude, offset, start_time)
    end

    eqs = [
        output.u ~ equation,
    ]

    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

# TODO:
# - Exponentials    Generate a rising and falling exponential signal
# - Pulse   Generate pulse signal of type Real
# - SawTooth    Generate saw tooth signal
# - Trapezoid   Generate trapezoidal signal of type Real
