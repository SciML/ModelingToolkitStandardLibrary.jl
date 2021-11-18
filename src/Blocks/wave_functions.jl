# Define and register smooth functions
using ModelingToolkit, Symbolics
_constant(K) = K
_cos_wave(x, f, A, st, ϕ) = A*cos(2*π*f*(x-st) + ϕ)
_damped_sine_wave(x, f, A, st, ϕ, d) = exp((st-x)*d)*A*sin(2*π*f*(x-st) + ϕ)
_ramp(x, δ, st, et, h) = h/(et-st)*(_xH(x, δ, st) - _xH(x, δ, et))
_sine(x, f, A, st, ϕ) = A*sin(2*π*f*(x - st) + ϕ)
_square_wave(x, δ, f, A, st) = A*2atan(sin(2π*(x-st)*f)/δ)/π
_step(x, δ, h, a) = h*(atan((x-a)/δ)/π + 0.5)
_triangular_wave(x, δ, f, A, st) = A*(1-2acos((1 - δ)sin(2π*(x-st)*f))/π)
_xH(x, δ, tₒ) = 0.5*(x-tₒ)*(1+((x-tₒ)/sqrt((x-tₒ)^2+δ^2)))


@register _constant(K)
@register _cos_wave(x, f, A, st, ϕ)
@register _damped_sine_wave(x, f, A, st, ϕ, damping)
@register _ramp(x, δ, st, et, h)
@register _sine(x, f, A, st, ϕ)
@register _square_wave(x, δ, f, A, st)
@register _step(x, δ, h, a)
@register _triangular_wave(x, δ, f, A, st)

function ConstantFunction(;name, K=1.0)
    val = K
    @info val
    @parameters K
    @variables y(t) [output=true]

    eqs = [
        y ~ val
    ]

    ODESystem(eqs, t, [y], [K], defaults=Dict(K => val), name=name)
end

function SmoothCosineFunction(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    o, A, f, st, ϕ = offset, amplitude, frequency, starttime, phase
    δ = 0.00001

    @parameters offset amplitude frequency starttime phase
    @variables y(t) [output=true]

    eqs = [
           y ~ _cos_wave(t, f, A, st, ϕ) * _step(t, δ, 1.0, st) + offset
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase), (o, A, f, st, ϕ)))
    ODESystem(eqs, t, [y], [offset, amplitude, frequency, starttime, phase], defaults=defaults, name=name)
end

function SmoothDampedSineFunction(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0, damping_coef=0.0)
    o, A, f, st, ϕ, d = offset, amplitude, frequency, starttime, phase, damping_coef
    δ = 0.0001

    @parameters offset amplitude frequency starttime phase damping_coef
    @variables y(t) [output=true]

    eqs = [
           y ~ _step(t, δ, o, 0.0) + _damped_sine_wave(t, f, A, st, ϕ, d) * _step(t, δ, 1.0, st)
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase, damping_coef), (o, A, f, st, ϕ, d)))
    ODESystem(eqs, t, [y], [offset, amplitude, frequency, starttime, phase, damping_coef], defaults=defaults, name=name)
end

function SmoothRampFunction(;name, offset=0.0, starttime=0.0, endtime=1.0, height=1.0)
    o, st, et, h = offset, starttime, endtime, height
    δ = 0.0001

    @parameters offset starttime endtime height
    @variables y(t) [output=true]

    eqs = [
           y ~ offset + _ramp(t, δ, st, et, h)
          ]
    defaults = Dict(zip((offset, starttime, endtime, height), (o, st, et, h)))
    ODESystem(eqs, t, [y], [offset, starttime, endtime, height], defaults=defaults, name=name)
end

function SmoothSineFunction(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
    o, A, f, st, ϕ = offset, amplitude, frequency, starttime, phase

    @parameters offset amplitude frequency starttime phase
    @variables y(t) [output=true]

    eqs = [
           y ~ offset + (t > st) * _sine(x, f, A, st, ϕ)
          ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime, phase), (o, A, f, st, ϕ)))
    ODESystem(eqs, t, [y], [offset, amplitude, frequency, starttime, phase], defaults=defaults, name=name)
end

function SmoothSquareFunction(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    o, A, f, st  = offset, amplitude, frequency, starttime
    δ = 0.0001

    @parameters offset amplitude frequency starttime
    @variables y(t) [output=true]

    eqs = [
        y ~ o + _square_wave(t, δ, f, A, st) * (t > st)
        ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime), (o, A, f, st)))
    ODESystem(eqs, t, [y], [offset, amplitude, frequency, starttime], defaults=defaults, name=name)
end

function SmoothStepFunction(;name, offset=0.0, starttime=0.0, height=1.0)
    o, st, h = offset, starttime, height
    δ = 0.0001

    @parameters offset starttime height
    @variables y(t)

    eqs = [
           y ~ offset + _step(t, δ, h, st)
          ]
    defaults = Dict(zip((offset, starttime, height), (o, st, h)))
    ODESystem(eqs, t, [y], [offset, starttime, height], defaults=defaults, name=name)
end

function SmoothTriangularFunction(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    o, A, f, st  = offset, amplitude, frequency, starttime
    δ = 0.0001

    @parameters offset amplitude frequency starttime
    @variables y(t) [output=true]
    
    eqs = [
        y ~ offset + (t>st) * _triangular_wave(t, δ, f, A, st)
    ]
    defaults = Dict(zip((offset, amplitude, frequency, starttime), (o, A, f, st)))
    ODESystem(eqs, t, [y], [offset, amplitude, frequency, starttime], defaults=defaults, name=name)
end
