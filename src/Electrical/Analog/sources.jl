# Voltage sources
function VoltageSource(;name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [v], [], systems=[p, n], name=name)
end

# Current Sources
function CurrentSource(;name)

    @named p = Pin()
    @named n = Pin()
    @variables i(t) v(t)

    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    ODESystem(eqs, t, [i, v], [], systems=[p, n], name=name)
end