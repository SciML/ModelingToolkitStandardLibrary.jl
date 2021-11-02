function Not(; name)
    @named x = DigitalPin()
    @named y = DigitalPin()

    eqs = [
        x.i ~ y.i
        y.val ~ _not(x.val)
    ]
    ODESystem(eqs, t, [], [], systems=[x, y], name=name)
end

function And(; name, N=2)
    x = map(1:N) do i
        DigitalPin(name=Symbol(:x, i))
    end
    @named y = DigitalPin()

    vals = [k.val for k in x]
    eqs = [
        y.val ~ _and(vals...)
        y.i ~ sum(k -> k.i, x)
    ]
    ODESystem(eqs, t, [], [], systems=[x..., y], name=name)
end

function Nand(; name, N=2)
    x = map(1:N) do i
        DigitalPin(name=Symbol(:x, i))
    end
    @named y = DigitalPin()

    vlist = [k.val for k in x]
    eqs = [
        y.val ~ _not(_and(vlist...))
        y.i ~ sum(k -> k.i, x) 
    ]
    ODESystem(eqs, t, [], [], systems=[x..., y], name=name)
end

function Or(; name, N=2)
    x = map(1:N) do i
        DigitalPin(name=Symbol(:x, i))
    end
    @named y = DigitalPin()

    vals = [k.val for k in x]
    eqs = [
        y.val ~ _or(vals...)
        y.i ~ sum(k -> k.i, x)
    ]
    ODESystem(eqs, t, [], [], systems=[x..., y], name=name)
end

function Nor(; name, N=2)
    x = map(1:N) do i
        DigitalPin(name=Symbol(:x, i))
    end
    @named y = DigitalPin()

    vlist = [k.val for k in x]
    eqs = [
        y.val ~ _not(_or(vlist...))
        y.i ~ sum(k -> k.i, x) 
    ]
    ODESystem(eqs, t, [], [], systems=[x..., y], name=name)
end

function Xor(; name, N=2)
    x = map(1:N) do i
        DigitalPin(name=Symbol(:x, i))
    end
    @named y = DigitalPin()

    vals = [k.val for k in x]
    eqs = [
        y.val ~ _xor(vals...)
        y.i ~ sum(k -> k.i, x)
    ]
    ODESystem(eqs, t, [], [], systems=[x..., y], name=name)
end

function Xnor(; name, N=2)
    x = map(1:N) do i
        DigitalPin(name=Symbol(:x, i))
    end
    @named y = DigitalPin()

    vlist = [k.val for k in x]
    eqs = [
        y.val ~ _not(_xor(vlist...))
        y.i ~ sum(k -> k.i, x) 
    ]
    ODESystem(eqs, t, [], [], systems=[x..., y], name=name)
end
