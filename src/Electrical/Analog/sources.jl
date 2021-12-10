# Define and register smooth functions
_cos_wave(x, f, A, st, ϕ) = A*cos(2*π*f*(x-st) + ϕ)
_damped_sine_wave(x, f, A, st, ϕ, d) = exp((st-x)*d)*A*sin(2*π*f*(x-st) + ϕ)
_ramp(x, δ, st, et, h) = h/(et-st)*(_xH(x, δ, st) - _xH(x, δ, et))
_square_wave(x, δ, f, A, st) = A*2atan(sin(2π*(x-st)*f)/δ)/π
_step(x, δ, h, a) = h*(atan((x-a)/δ)/π + 0.5)
_triangular_wave(x, δ, f, A, st) = A*(1-2acos((1 - δ)sin(2π*(x-st)*f))/π)
_xH(x, δ, tₒ) = 0.5*(x-tₒ)*(1+((x-tₒ)/sqrt((x-tₒ)^2+δ^2)))

@register _cos_wave(x, f, A, st, ϕ)
@register _damped_sine_wave(x, f, A, st, ϕ, damping)
@register _ramp(x, δ, st, et, h)
@register _square_wave(x, δ, f, A, st)
@register _step(x, δ, h, a)
@register _triangular_wave(x, δ, f, A, st)

# Voltage sources
"""
```julia
function ConstantVoltage(;name, V=1.0)
```

The source for an ideal constant voltage.

# Observables
- `V`
  The constant voltage across the terminals of this source

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v` and is always constant

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin
"""
function ConstantVoltage(;name, V=1.0)
    val = V

    @named p = Pin()
    @named n = Pin()
    @parameters V
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ V
          ]
    ODESystem(eqs, t, [v], [V], systems=[p, n], defaults=Dict(V => val), name=name)
end

"""
```julia
function CosineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
```

A source in which the voltage across its terminals is a cosine function of time.

# Observables
- `offset`
  A constant offset added to the voltage output
- `amplitude`
  The amplitude of the cosine function
- `frequency`
  The frequency of the cosine function
- `starttime`
  The time at which the source starts functioning. Before this time, the voltage across
  its terminals is 0.
- `phase`
  The phase offset of the cosine function

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function CosineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    o, A, f, st, ϕ = offset, amplitude, frequency, starttime, phase
    δ = 0.00001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime phase
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           v ~ _cos_wave(t, f, A, st, ϕ) * _step(t, δ, 1.0, st) + offset
           0 ~ p.i + n.i
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase), (o, A, f, st, ϕ)))
    ODESystem(eqs, t, [v], [offset, amplitude, frequency, starttime, phase], systems=[p, n], defaults=defaults, name=name)
end

"""
```julia
function DampedSineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0, damping_coef=0.0)
```

A source in which the voltage across its terminals is a damped sine function of time.

# Observables
- `offset`
  A constant offset added to the voltage output
- `amplitude`
  The amplitude of the damped sine function
- `frequency`
  The frequency of the damped sine function
- `starttime`
  The time at which the source starts functioning. Before this time, the voltage across
  its terminals is 0.
- `phase`
  The phase offset of the damped sine function
- `damping_coef`
  Damping coefficient of the damped sine function

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function DampedSineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0, damping_coef=0.0)
    o, A, f, st, ϕ, d = offset, amplitude, frequency, starttime, phase, damping_coef
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime phase damping_coef
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           v ~ _step(t, δ, o, 0.0) + _damped_sine_wave(t, f, A, st, ϕ, d) * _step(t, δ, 1.0, st)
           0 ~ p.i + n.i
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase, damping_coef), (o, A, f, st, ϕ, d)))
    ODESystem(eqs, t, [v], [offset, amplitude, frequency, starttime, phase, damping_coef], systems=[p, n], defaults=defaults, name=name)
end

"""
```julia
function RampVoltage(;name, offset=0.0, starttime=0.0, endtime=1.0, height=1.0)
```

A source in which the voltage across grows linearly from `offset` to `offset+height` over
the time interval `[starttime, endtime]`

# Observables
- `offset`
  A constant offset added to the voltage output
- `starttime`
  The time at which the voltage starts growing
- `endtime`
  The time at which the voltage finishes growing
- `height`
  The amount that the voltage grows in the interval `[starttime, endtime]`

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function RampVoltage(;name, offset=0.0, starttime=0.0, endtime=1.0, height=1.0)
    o, st, et, h = offset, starttime, endtime, height
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset starttime endtime height
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           v ~ offset + _ramp(t, δ, st, et, h)
           0 ~ p.i + n.i
          ]
    defaults = Dict(zip((offset, starttime, endtime, height), (o, st, et, h)))
    ODESystem(eqs, t, [v], [offset, starttime, endtime, height], systems=[p, n], defaults=defaults, name=name)
end

"""
```julia
function SineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
```

A source in which the voltage across its terminals is a sine function of time.

# Observables
- `offset`
  A constant offset added to the voltage output
- `amplitude`
  The amplitude of the sine function
- `frequency`
  The frequency of the sine function
- `starttime`
  The time at which the source starts functioning. Before this time, the voltage across
  its terminals is 0.
- `phase`
  The phase offset of the sine function

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function SineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    o, A, f, st, ϕ = offset, amplitude, frequency, starttime, phase

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime phase
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           v ~ offset + (t > st) * (A*sin(2*π*f*(t - st) + ϕ))
           0 ~ p.i + n.i
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase), (o, A, f, st, ϕ)))
    ODESystem(eqs, t, [v], [offset, amplitude, frequency, starttime, phase], systems=[p, n], defaults=defaults, name=name)
end

"""
```julia
function SquareVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
```

A source in which the voltage across its terminals is a square function of time.

# Observables
- `offset`
  A constant offset added to the voltage output
- `amplitude`
  The amplitude of the square wave function
- `frequency`
  The frequency of the square wave function
- `starttime`
  The time at which the source starts functioning. Before this time, the voltage across
  its terminals is 0.

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function SquareVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    o, A, f, st  = offset, amplitude, frequency, starttime
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime
    @variables v(t)
    
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        v ~ o + _square_wave(t, δ, f, A, st) * (t > st)
        ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime), (o, A, f, st)))
    ODESystem(eqs, t, [v], [offset, amplitude, frequency, starttime], systems=[p, n], defaults=defaults, name=name)
end

"""
```julia
function StepVoltage(;name, offset=0.0, starttime=0.0, height=1.0)
```

A source in which the voltage across its terminals increases from `offset` to `offset+height` at
`starttime`

# Observables
- `offset`
  A constant offset added to the voltage output
- `starttime`
  The time at which the source starts functioning, and the voltage jumps
- `height`
  Magnitude of increase in voltage

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function StepVoltage(;name, offset=0.0, starttime=0.0, height=1.0)
    o, st, h = offset, starttime, height
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset starttime height
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           v ~ offset + _step(t, δ, h, st)
           0 ~ p.i + n.i
          ]
    defaults = Dict(zip((offset, starttime, height), (o, st, h)))
    ODESystem(eqs, t, [v], [offset, starttime, height], systems=[p, n], defaults=defaults, name=name)
end

"""
```julia
function TriangularVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
```

A source in which the voltage across its terminals is a triangular function of time.

# Observables
- `offset`
  A constant offset added to the voltage output
- `amplitude`
  Amplitude of the triangular wave function
- `frequency`
  Frequency of the triangular wave function
- `starttime`
  The time at which the source starts functioning. Before this, the output of the source is 0

# States
- `v(t)`
  The voltage across this source, given by `p.v - n.v`

# Connectors
- `p`
  Positive port
- `n`
  Negative port
"""
function TriangularVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    o, A, f, st  = offset, amplitude, frequency, starttime
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime
    @variables v(t)
    
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        v ~ offset + (t>st) * _triangular_wave(t, δ, f, A, st)
    ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime), (o, A, f, st)))
    ODESystem(eqs, t, [v], [offset, amplitude, frequency, starttime], systems=[p, n], defaults=defaults, name=name)
end

# Current Sources
function ConstantCurrent(;name, I=1.0)
    val = I

    @named p = Pin()
    @named n = Pin()
    @parameters I
    @variables i(t)

    eqs = [
           0 ~ p.i + n.i
           i ~ p.i
           i ~ I
          ]
    ODESystem(eqs, t, [i], [I], systems=[p, n], defaults=Dict(I => val), name=name)
end

function CosineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    o, A, f, st, ϕ = offset, amplitude, frequency, starttime, phase
    δ = 0.00001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime phase
    @variables i(t)

    eqs = [
           i ~ _cos_wave(t, f, A, st, ϕ) * _step(t, δ, 1.0, st) + offset
           0 ~ p.i + n.i
           i ~ p.i
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase), (o, A, f, st, ϕ)))
    ODESystem(eqs, t, [i], [offset, amplitude, frequency, starttime, phase], systems=[p, n], defaults=defaults, name=name)
end

function DampedSineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0, damping_coef=1.0)
    o, A, f, st, ϕ, d = offset, amplitude, frequency, starttime, phase, damping_coef
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime phase damping_coef
    @variables i(t)

    eqs = [
           i ~ _step(t, δ, o, 0.0) + _damped_sine_wave(t, f, A, st, ϕ, d) * _step(t, δ, 1.0, st)
           0 ~ p.i + n.i
           i ~ p.i
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase, damping_coef), (o, A, f, st, ϕ, d)))
    ODESystem(eqs, t, [i], [offset, amplitude, frequency, starttime, phase, damping_coef], systems=[p, n], defaults=defaults, name=name)
end

function RampCurrent(;name, offset=0.0, starttime=0.0, endtime=1.0, height=1.0)
    o, st, et, h = offset, starttime, endtime, height
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset starttime endtime height
    @variables i(t)

    eqs = [
           i ~ _step(t, δ, o, 0.0) + _ramp(t, δ, st, et, h)
           0 ~ p.i + n.i
           i ~ p.i
          ]
    defaults = Dict(zip((offset, starttime, endtime, height), (o, st, et, h)))
    ODESystem(eqs, t, [i], [offset, starttime, endtime, height], systems=[p, n], defaults=defaults, name=name)
end

function SineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    o, A, f, st, ϕ = offset, amplitude, frequency, starttime, phase

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime phase
    @variables i(t)

    eqs = [
           i ~ offset + (t > st) * (A*sin(2*π*f*(t - st) + ϕ))
           0 ~ p.i + n.i
           i ~ p.i
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase), (o, A, f, st, ϕ)))
    ODESystem(eqs, t, [i], [offset, amplitude, frequency, starttime, phase], systems=[p, n], defaults=defaults, name=name)
end

function SquareCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    o, A, f, st  = offset, amplitude, frequency, starttime
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime
    @variables i(t)
    
    eqs = [
        0 ~ p.i + n.i
        i ~ o + _square_wave(t, δ, f, A, st) * (t > st)
        i ~ p.i
        ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime), (o, A, f, st)))
    ODESystem(eqs, t, [i], [offset, amplitude, frequency, starttime], systems=[p, n], defaults=defaults, name=name)
end

function StepCurrent(;name, offset=0.0, starttime=0.0, height=1.0)
    o, st, h = offset, starttime, height
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset starttime height
    @variables i(t)

    eqs = [
           i ~ offset + _step(t, δ, h, st)
           0 ~ p.i + n.i
           i ~ p.i
          ]
    defaults = Dict(zip((offset, starttime, height), (o, st, h)))
    ODESystem(eqs, t, [i], [offset, starttime, height], systems=[p, n], defaults=defaults, name=name)
end

function TriangularCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    o, A, f, st  = offset, amplitude, frequency, starttime
    δ = 0.0001

    @named p = Pin()
    @named n = Pin()
    @parameters offset amplitude frequency starttime
    @variables i(t)
    
    eqs = [
        0 ~ p.i + n.i
        i ~ offset + _step(t, δ, 1, st) * _triangular_wave(t, δ, f, A, st)
        i ~ p.i
    ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime), (o, A, f, st)))
    ODESystem(eqs, t, [i], [offset, amplitude, frequency, starttime], systems=[p, n], defaults=defaults, name=name)
end
