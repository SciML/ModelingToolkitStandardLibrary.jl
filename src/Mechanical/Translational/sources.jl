#=
function SineForce(; name, amp, freq)
    @named flange = MechanicalPort()

    pars = @parameters begin
        amp = amp
        freq = freq
    end
    vars = []
    eqs = [
        flange.f ~ amp * sin(2 * π * t * freq),
    ]
    compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.v => 0]), flange)
end
=#
