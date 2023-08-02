using DiffEqBase
import ChainRulesCore

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
    Constant(; name, k = 0.0)

Generate constant signal.

# Parameters:

  - `k`: Constant output value

# Connectors:

  - `output`
"""
@mtkmodel Constant begin
    @components begin
        output = RealOutput()
    end
    @parameters begin
        k = 0.0, [description = "Constant output value of block"]
    end
    @equations begin
        output.u ~ k
    end
end

"""
    TimeVaryingFunction(f; name)

Outputs ``f(t)``.

The input variable `t` can be changed by passing a different variable as the keyword argument `t`.

# Connectors:
- `output`
"""
@mtkmodel TimeVaryingFunction begin
    @parameters begin
        f
    end
    @components begin
        output = RealOutput()
    end
    @equations begin
        output.u ~ first(getdefault(f))(t)
    end
end
TimeVaryingFunction.f(f; name) = TimeVaryingFunction.f(; f = [f], name)

"""
    Sine(; name, frequency, amplitude = 1, phase = 0, offset = 0, start_time = 0,
    smooth = false)

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
    Cosine(; name, frequency, amplitude = 1, phase = 0, offset = 0, start_time = 0,
    smooth = false)

Cosine signal.

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
    ContinuousClock(; name, offset = 0, start_time = 0)

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
Ramp(; name, height = 1, duration = 1, offset = 0, start_time = 0, smooth = false)

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
    Square(; name, frequency = 1.0, amplitude = 1.0, offset = 0.0, start_time = 0.0,
    smooth = false)
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
    ExpSine(; name, frequency, amplitude = 1, damping = 0.1, phase = 0, offset = 0, start_time = 0, smooth = false)

Exponentially damped sine signal.

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
    Triangular(; name, amplitude = 1.0, frequency = 1.0, offset = 0.0,
    start_time = 0.0, smooth = false)

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

function linear_interpolation(x1::T, x2::T, t1::T, t2::T, t) where {T <: Real}
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
    circular_buffer::Bool
end

Parameter(data::Vector{T}, ref::T) where {T <: Real} = Parameter(data, ref, true)
Parameter(x::Parameter) = x
function Parameter(x::T; tofloat = true) where {T <: Real}
    if tofloat
        x = float(x)
        P = typeof(x)
    else
        P = T
    end

    return Parameter(P[], x)
end

function Base.isequal(x::Parameter, y::Parameter)
    b0 = length(x.data) == length(y.data)
    if b0
        b1 = all(x.data .== y.data)
        b2 = x.ref == y.ref
        return b1 & b2
    else
        return false
    end
end

Base.:*(x::Number, y::Parameter) = x * y.ref
Base.:*(y::Parameter, x::Number) = Base.:*(x, y)
Base.:*(x::Parameter, y::Parameter) = x.ref * y.ref

Base.:/(x::Number, y::Parameter) = x / y.ref
Base.:/(y::Parameter, x::Number) = y.ref / x
Base.:/(x::Parameter, y::Parameter) = x.ref / y.ref

Base.:+(x::Number, y::Parameter) = x + y.ref
Base.:+(y::Parameter, x::Number) = Base.:+(x, y)
Base.:+(x::Parameter, y::Parameter) = x.ref + y.ref

Base.:-(y::Parameter) = -y.ref
Base.:-(x::Number, y::Parameter) = x - y.ref
Base.:-(y::Parameter, x::Number) = y.ref - x
Base.:-(x::Parameter, y::Parameter) = x.ref - y.ref

Base.:^(x::Number, y::Parameter) = Base.:^(x, y.ref)
Base.:^(y::Parameter, x::Number) = Base.:^(y.ref, x)
Base.:^(x::Parameter, y::Parameter) = Base.:^(x.ref, y.ref)

Base.isless(x::Parameter, y::Number) = Base.isless(x.ref, y)
Base.isless(y::Number, x::Parameter) = Base.isless(y, x.ref)

Base.copy(x::Parameter{T}) where {T} = Parameter{T}(copy(x.data), x.ref)

ifelse(c::Bool, x::Parameter, y::Parameter) = ifelse(c, x.ref, y.ref)
ifelse(c::Bool, x::Parameter, y::Number) = ifelse(c, x.ref, y)
ifelse(c::Bool, x::Number, y::Parameter) = ifelse(c, x, y.ref)

Base.max(x::Number, y::Parameter) = max(x, y.ref)
Base.max(x::Parameter, y::Number) = max(x.ref, y)
Base.max(x::Parameter, y::Parameter) = max(x.ref, y.ref)

Base.min(x::Number, y::Parameter) = min(x, y.ref)
Base.min(x::Parameter, y::Number) = min(x.ref, y)
Base.min(x::Parameter, y::Parameter) = min(x.ref, y.ref)

function Base.show(io::IO, m::MIME"text/plain", p::Parameter)
    if !isempty(p.data)
        print(io, p.data)
    else
        print(io, p.ref)
    end
end

function get_sampled_data(t, memory::Parameter{T}) where {T}
    if t < 0
        t = zero(t)
    end

    if isempty(memory.data)
        if T <: AbstractFloat
            return T(NaN)
        else
            return zero(T)
        end
    end

    i1 = floor(Int, t / memory.ref) + 1 #expensive
    i2 = i1 + 1

    t1 = (i1 - 1) * memory.ref
    x1 = @inbounds memory.data[i1]

    if t == t1
        return x1
    else
        n = length(memory.data)

        if memory.circular_buffer
            i1 = (i1 - 1) % n + 1
            i2 = (i2 - 1) % n + 1
        else
            if i2 > n
                i2 = n
                i1 = i2 - 1
            end
        end

        t2 = (i2 - 1) * memory.ref
        x2 = @inbounds memory.data[i2]
        return linear_interpolation(x1, x2, t1, t2, t)
    end
end

get_sample_time(memory::Parameter) = memory.ref
Symbolics.@register_symbolic get_sample_time(memory)

Symbolics.@register_symbolic get_sampled_data(t, memory)

function first_order_backwards_difference(t, memory)
    Δt = get_sample_time(memory)
    x1 = get_sampled_data(t, memory)
    x0 = get_sampled_data(t - Δt, memory)

    return (x1 - x0) / Δt
end

function Symbolics.derivative(::typeof(get_sampled_data), args::NTuple{2, Any}, ::Val{1})
    t = @inbounds args[1]
    memory = @inbounds args[2]
    first_order_backwards_difference(t, memory)
end
function ChainRulesCore.frule((_, ẋ, _), ::typeof(get_sampled_data), t, memory)
    first_order_backwards_difference(t, memory) * ẋ
end

"""
    SampledData(; name, buffer)

data input component.

# Parameters:
  - `buffer`: a `Parameter` type which holds the data and sample time

# Connectors:
  - `output`
"""
@component function SampledData(; name, buffer)
    pars = @parameters begin
        buffer = buffer
    end
    vars = []
    systems = @named begin
        output = RealOutput()
    end
    eqs = [
        output.u ~ get_sampled_data(t, buffer),
    ]
    return ODESystem(eqs, t, vars, pars; name, systems,
        defaults = [output.u => get_sampled_data(0.0, buffer)])
end
@deprecate Input SampledData

function SampledData(T::Type, circular_buffer = true; name)
    SampledData(T[], zero(T), circular_buffer; name)
end
function SampledData(dt::T, circular_buffer = true) where {T <: Real}
    SampledData(T[], dt, circular_buffer; name)
end
function SampledData(data::Vector{T}, dt::T, circular_buffer = true; name) where {T <: Real}
    SampledData(; name, buffer = Parameter(data, dt, circular_buffer))
end

Base.convert(::Type{T}, x::Parameter{T}) where {T <: Real} = x.ref
function Base.convert(::Type{<:Parameter{T}}, x::Number) where {T <: Real}
    Parameter{T}(T[], x, true)
end

# Beta Code for potential AE Hack ----------------------
function set_sampled_data!(memory::Parameter{T}, t, x, Δt::Parameter{T}) where {T}
    if t < 0
        t = zero(t)
    end

    if t == zero(t)
        empty!(memory.data)
    end

    n = length(memory.data)
    i = round(Int, t / Δt) + 1 #expensive
    if i == n + 1
        push!(memory.data, DiffEqBase.value(x))
    elseif i <= n
        @inbounds memory.data[i] = DiffEqBase.value(x)
    else
        error("Memory buffer skipped a step: n=$n, i=$i")
    end

    # memory.ref = Δt

    return x
end
Symbolics.@register_symbolic set_sampled_data!(memory, t, x, Δt)

function Symbolics.derivative(::typeof(set_sampled_data!), args::NTuple{4, Any}, ::Val{2})
    memory = @inbounds args[1]
    t = @inbounds args[2]
    x = @inbounds args[3]
    Δt = @inbounds args[4]
    first_order_backwards_difference(t, x, Δt, memory)
end
Symbolics.derivative(::typeof(set_sampled_data!), args::NTuple{4, Any}, ::Val{3}) = 1 #set_sampled_data returns x, therefore d/dx (x) = 1
function ChainRulesCore.frule((_, _, ṫ, ẋ, _),
    ::typeof(set_sampled_data!),
    memory,
    t,
    x,
    Δt)
    first_order_backwards_difference(t, x, Δt, memory) * ṫ + ẋ
end

function first_order_backwards_difference(t, x, Δt, memory)
    x1 = set_sampled_data!(memory, t, x, Δt)
    x0 = get_sampled_data(t - Δt, memory)

    return (x1 - x0) / Δt
end
