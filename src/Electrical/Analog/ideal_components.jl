"""
    Ground(;name)

Ground of an electrical circuit. The potential at the ground node is zero. Every circuit must have one ground node.
"""
function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], []; systems=[g], name=name)
end

"""
    Resistor(;name, R=1.0)

Ideal linear electrical resistor.

# Parameters: 
- `R`: [Ω] Resistance
"""
function Resistor(;name, R=1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters R=R
    eqs = [
        v ~ i * R
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    Capacitor(;name, C=1.0, v_start=0.0) 

Ideal linear electrical capacitor.

# Parameters:
- `C`: [F] Capacity
- `v_start`: [V] Initial voltage of capacitor
"""
function Capacitor(;name, C=1.0, v_start=0.0) 
    @named oneport = OnePort(;v_start=v_start)
    @unpack v, i = oneport
    pars = @parameters C=C
    eqs = [
        D(v) ~ i / C
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    Inductor(;name, L=1.0e-6, i_start=0.0)

Ideal linear electrical inductor.

# Parameters:
- `L`: [H] Inductance
- `i_start`: [A] Initial current through inductor
"""
function Inductor(;name, L=1.0e-6, i_start=0.0)
    @named oneport = OnePort(;i_start=i_start)
    @unpack v, i = oneport
    pars = @parameters L=L
    eqs = [
        D(i) ~ 1 / L * v
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    IdealOpAmp(;name)

Ideal operational amplifier (norator-nullator pair).
The ideal OpAmp is a two-port. The left port is fixed to v1=0 and i1=0 (nullator). 
At the right port both any voltage v2 and any current i2 are possible (norator).
"""
function IdealOpAmp(;name)
    @named p1 = Pin()
    @named p2 = Pin()
    @named n1 = Pin()
    @named n2 = Pin()
    sts = @variables begin
        v1(t) 
        v2(t) 
        i1(t) 
        i2(t)
    end
    eqs = [
        v1 ~ p1.v - n1.v
        v2 ~ p2.v - n2.v
        0 ~ p1.i + n1.i
        0 ~ p2.i + n2.i
        i1 ~ p1.i
        i2 ~ p2.i
        v1 ~ 0
        i1 ~ 0
    ]
    ODESystem(eqs, t, sts, [], systems=[p1, p2, n1, n2], name=name)
end
