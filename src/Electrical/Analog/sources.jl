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
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency
        starttime=starttime
        phase=phase 
        damping_coef=damping_coef
    end
    eqs = [
        v ~ _damped_sine_wave(t, frequency, amplitude, starttime, phase, damping_coef) * _step(t, δ, 1.0, starttime)
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
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

function SineVoltage(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
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
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        starttime=starttime 
    end
    eqs = [
        v ~ _step(t, δ, height, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function TriangularVoltage(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        starttime=starttime 
    end
    eqs = [
        v ~ _triangular_wave(t, δ, frequency, amplitude, starttime) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

# Current Sources
function ConstantCurrent(;name, 
    I = 1.0, # [A]
    )   
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters A=A
    eqs = [
        i ~ I
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function CosineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
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
        i ~ _cos_wave(t, frequency, amplitude, starttime, phase) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function DampedSineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0, damping_coef=0.0)
    δ = 0.00001
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency
        starttime=starttime
        phase=phase 
        damping_coef=damping_coef
    end
    eqs = [
        i ~ _damped_sine_wave(t, frequency, amplitude, starttime, phase, damping_coef) * _step(t, δ, 1.0, starttime)
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function RampCurrent(;name, offset=0.0, starttime=0.0, endtime=1.0, height=1.0)
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
        i ~ _ramp(t, δ, 1.0, starttime, height) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function SineCurrent(;name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0, phase=0.0)
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
        i ~ _sin_wave(t, frequency, amplitude, starttime, phase) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function SquareCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        starttime=starttime 
    end
    eqs = [
        i ~ _square_wave(t, frequency, amplitude, starttime, phase) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function StepCurrent(;name, offset=0.0, starttime=0.0, height=1.0)
    δ = 0.0001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        height=height 
        starttime=starttime 
    end
    eqs = [
        i ~ _step(t, δ, height, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

function TriangularCurrent(; name, offset=0.0, amplitude=1.0, frequency=1.0, starttime=0.0)
    δ = 0.00001

    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters begin
        offset=offset
        amplitude=amplitude 
        frequency=frequency 
        starttime=starttime 
    end
    eqs = [
        i ~ _triangular_wave(t, δ, frequency, amplitude, starttime) * _step(t, δ, 1.0, starttime) + offset
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end
