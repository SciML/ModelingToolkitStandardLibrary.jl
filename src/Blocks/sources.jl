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
    pars = @parameters offset=offset startTime=starttime amplitude=amplitude frequency=frequency phase=phase
    eqs = [
        output.u ~ offset + ifelse(t < startTime, 0, amplitude* sin(2*pi*frequency*(t - startTime) + phase))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [output])
end

# TODO:
# - Clock   Generate actual time signal
# - Step    Generate step signal of type Real
# - Ramp    Generate ramp signal
# - ExpSine Generate exponentially damped sine signal
# - Exponentials    Generate a rising and falling exponential signal
# - Pulse   Generate pulse signal of type Real
# - SawTooth    Generate saw tooth signal
# - Trapezoid   Generate trapezoidal signal of type Real