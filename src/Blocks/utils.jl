@connector function RealInput(;name)
    sts = @variables u(t)=0
    ODESystem(Equation[], t, sts, []; name=name)
end

@connector function RealOutput(;name)
    sts = @variables u(t)=0
    ODESystem(Equation[], t, sts, []; name=name)
end

function SISO(;name)
    @named u = RealInput()
    @named y = RealOutput()
    return ODESystem(; name, systems=[u, y])
end