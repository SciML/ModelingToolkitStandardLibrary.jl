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

# Voltage sources
function ConstantVoltage(;name, 
    V = 1.0, # [V]
    )   
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters V=V
    eqs = [
        v ~ V
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function CosineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        starttime=starttime 
        phase=phase
    end
    eqs = [
        v ~ _cos_wave(t, frequency, amplitude, starttime, phase) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

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

function RampVoltage(;name, offset=0.0, starttime=0.0, endtime=1.0, height=1.0)
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        starttime=starttime 
        endtime=endtime
    end
    eqs = [
        v ~ _ramp(t, δ, 1.0, starttime, height) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function CosineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        starttime=starttime 
        phase=phase
    end
    eqs = [
        v ~ _sin_wave(t, frequency, amplitude, starttime, phase) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function SquareVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        starttime=starttime 
        endtime=endtime
    end
    eqs = [
        v ~ _square_wave(t, frequency, amplitude, starttime, phase) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

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
