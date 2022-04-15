"""
Generate constant signal.
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
"""
function SinSource(;name, 
    frequency,# [Hz] Frequency of sine wave 
    amplitude=1, # Amplitude of sine wave
    phase=0, # [rad] Phase of sine wave 
    offset=0, # Offset of output signal
    starttime=0)

    @named output = RealOutput()
    pars = @parameters offset=offset startime=starttime amplitude=amplitude frequency=frequency phase=phase
    eqs = [
        output.u ~ offset + ifelse(t < startime, 0, amplitude* sin(2*pi*frequency*(t - startime) + phase))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate clock signal.
"""
function ClockSource(;name, 
    offset=0, # Offset of output signal
    starttime=0)

    @named output = RealOutput()
    pars = @parameters offset=offset starttime=starttime
    eqs = [
        output.u ~ offset + ifelse(t < starttime, 0, t - starttime)
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate ramp signal.
"""
function RampSource(;name, 
    offset=0, # Offset of output signal
    height=1,
    duration, 
    starttime=0)

    @named output = RealOutput()
    pars = @parameters offset=offset starttime=starttime height=height duration=duration
    eqs = [
        output.u ~ offset + ifelse(t < starttime, 0, 
            ifelse(t < (starttime + duration), (t - starttime) * height / duration, height))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

"""
Generate step signal.
"""
function StepSource(;name, 
    offset=0, # Offset of output signal
    height=1,
    starttime=0)

    @named output = RealOutput()
    pars = @parameters offset=offset starttime=starttime height=height
    eqs = [
        output.u ~ offset + ifelse(t < starttime, 0, height)
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

# TODO:
# - ExpSine Generate exponentially damped sine signal
# - Exponentials    Generate a rising and falling exponential signal
# - Pulse   Generate pulse signal of type Real
# - SawTooth    Generate saw tooth signal
# - Trapezoid   Generate trapezoidal signal of type Real