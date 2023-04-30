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
@component function Constant(; name, k = 1)
    @named output = RealOutput()
    pars = @parameters k=k [description = "Constant output value of block $name"]
    eqs = [
        output.u ~ k,
    ]
    compose(ODESystem(eqs, t, [], pars; name = name), [output])
end

"""
    TimeVaryingFunction(f; t=t, name)

Outputs ``f(t)``.

The input variable `t` can be changed by passing a different variable as the keyword argument `t`.

# Connectors:
- `output`
"""
@component function TimeVaryingFunction(f; t = t, name)
    @named output = RealOutput()
    eqs = [
        output.u ~ f(t),
    ]
    compose(ODESystem(eqs, Blocks.t; name = name), [output])
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
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Sine(; name,
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
        smooth === true && (smooth = 1e-5)
        smooth_sin(t, smooth, frequency, amplitude, phase, offset, start_time)
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
             It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:
- `output`
"""

@component function Cosine(; name,
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
        smooth === true && (smooth = 1e-5)
        smooth_cos(t, smooth, frequency, amplitude, phase, offset, start_time)
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
@component function ContinuousClock(; name, offset = 0, start_time = 0)
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
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Ramp(; name,
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
        smooth === true && (smooth = 1e-5)
        smooth_ramp(t, smooth, height, duration, offset, start_time)
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
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Square(; name, frequency = 1.0, amplitude = 1.0,
                           offset = 0.0, start_time = 0.0, smooth = false)
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
        smooth === true && (smooth = 1e-5)
        smooth_square(t, smooth, frequency, amplitude, offset, start_time)
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
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `true`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Step(; name, height = 1, offset = 0, start_time = 0, duration = Inf,
                         smooth = 1e-5)
    @named output = RealOutput()
    duration_numeric = duration
    pars = @parameters offset=offset start_time=start_time height=height duration=duration
    equation = if smooth == false # use comparison in case smooth is a float
        offset + ifelse((start_time < t) & (t < start_time + duration), height, 0)
    else
        smooth === true && (smooth = 1e-5)
        if duration_numeric == Inf
            smooth_step(t, smooth, height, offset, start_time)
        else
            smooth_step(t, smooth, height, offset, start_time) -
            smooth_step(t, smooth, height, 0, start_time + duration)
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
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function ExpSine(; name,
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
        smooth === true && (smooth = 1e-5)
        smooth_damped_sin(t, smooth, frequency, amplitude, damping, phase, offset,
                          start_time)
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
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Triangular(; name, amplitude = 1.0, frequency = 1.0,
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
        smooth === true && (smooth = 1e-5)
        smooth_triangular(t, smooth, frequency, amplitude, offset, start_time)
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

function linear_interpolation(x1::T, x2::T, t1::T, t2::T, t::T) where {T <: Real}
    if t1 != t2
        slope = (x2 - x1) / (t2 - t1)
        intercept = x1 - slope * t1

        return slope * t + intercept
    else
        @assert x1==x2 "x1 ($x1) and x2 ($x2) should be equal if t1 == t2"

        return x2
    end
end

struct Parameter{T <: Real}
    data::Vector{T}
    ref::T
    n::Int
end

function Base.isequal(x::Parameter, y::Parameter) 
    b0 = x.n == y.n
    if b0
        b1 = all(x.data .== y.data)
        b2 = x.ref == y.ref
        return b1 & b2
    else
        return false
    end
end

Base.:*(x::Number, y::Parameter) = x*y.ref
Base.:*(y::Parameter, x::Number) = Base.:*(x,y)
Base.:*(x::Parameter, y::Parameter) = x.ref*y.ref

Base.:/(x::Number, y::Parameter) = x/y.ref
Base.:/(y::Parameter, x::Number) = y.ref/x
Base.:/(x::Parameter, y::Parameter) = x.ref/y.ref

Base.:+(x::Number, y::Parameter) = x + y.ref
Base.:+(y::Parameter, x::Number) = Base.:+(x,y)
Base.:+(x::Parameter, y::Parameter) = x.ref + y.ref

Base.:-(x::Number, y::Parameter) = x - y.ref
Base.:-(y::Parameter, x::Number) = y.ref - x
Base.:-(x::Parameter, y::Parameter) = x.ref - y.ref

Base.:^(x::Number, y::Parameter) = Base.power_by_squaring(x, y.ref)
Base.:^(y::Parameter, x::Number) = Base.power_by_squaring(y.ref, x)
Base.:^(x::Parameter, y::Parameter) = Base.power_by_squaring(x.ref, y.ref)

Base.isless(x::Parameter, y::Number) = Base.isless(x.ref, y)
Base.isless(y::Number, x::Parameter) = Base.isless(y, x.ref)




Base.copy(x::Parameter{T}) where {T} = Parameter{T}(copy(x.data), x.ref, x.n)

function Base.show(io::IO, m::MIME"text/plain", p::Parameter)
    if !isempty(p.data)
	    print(io, p.data)
    else
        print(io, p.ref)
    end
end

Parameter(x::Parameter) = x
function Parameter(x::T; tofloat=true) where T <: Real 
    if tofloat
        x = float(x)
        P = typeof(x)
    else
        P = T
    end

    return Parameter(P[], x, 0)
end
Parameter(x::Vector{T}, dt::T) where T <: Real = Parameter(x, dt, length(x))


function input(t, memory::Parameter)
    if t < 0
        t = zero(t)
    end

    i1 = floor(Int, t / memory.ref) + 1 #expensive
    i2 = i1 + 1

    if i2 > memory.n
        i2 = memory.n
        i1 = i2-1
    end

    t1 = i1 * memory.ref
    x1 = @inbounds getindex(memory.data, i1)

    if t == t1
        return x1
    else

        t2 = i2 * memory.ref
        x2 = @inbounds getindex(memory.data, i2)
        return linear_interpolation(x1, x2, t1, t2, t)
    end
    
end

get_sample_time(memory::Parameter) = memory.ref
Symbolics.@register_symbolic get_sample_time(memory)

Symbolics.@register_symbolic input(t, memory)

function first_order_backwards_difference(t, memory)
    Δt = get_sample_time(memory)
    x1 = input(t, memory)
    x0 = input(t - Δt, memory)

    return (x1 - x0) / Δt
end

function Symbolics.derivative(::typeof(input), args::NTuple{2, Any}, ::Val{1})
    first_order_backwards_difference(args[1], args[2])
end

Input(T::Type; name) = Input(T[], zero(T); name)
function Input(data::Vector{T}, dt::T; name) where {T <: Real}
    Input(; name, buffer = Parameter(data, dt))
end

"""
    Input(; name, buffer)

data input component.  

# Parameters:
  - `buffer`: a `Parameter` type which holds the data and sample time

# Connectors:
  - `output`
"""
@component function Input(; name, buffer)
    pars = @parameters begin buffer = buffer end
    vars = [] 
    systems = @named begin output = RealOutput() end
    eqs = [
        output.u ~ input(t, buffer)
    ]
    return ODESystem(eqs, t, vars, pars; name, systems, defaults=[output.u => 0.0]) #TODO: get initial value from buffer
end
