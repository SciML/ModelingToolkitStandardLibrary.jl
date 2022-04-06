@connector function RealInput(;name, nin=1)
    if nin == 1
        @variables u(t) = 0.0
    else
        @variables u[1:nin](t) = zeros(nin)
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name=name)
end

@connector function RealOutput(;name, nout=1)
    if nout == 1
        @variables u(t) = 0.0
    else
        @variables u[1:nout](t) = zeros(nout)
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name=name)
end

function SISO(;name)
    @named input = RealInput()
    @named output = RealOutput()
    @variables begin
        u(t)=0.0
        y(t)=0.0
    end
    eqs = [
        u ~ input.u
        y ~ output.u
    ]
    return ODESystem(eqs, t, [u, y], []; name, systems=[input, output])
end

function MIMO(;name, nin=1, nout=1)
    @named input = RealInput(nin=nin)
    @named output = RealOutput(nout=nout)
    @variables begin
        u[1:nin](t)=zeros(nin)
        y[1:nout](t)=zeros(nout)
    end
    eqs = [
        u ~ input.u
        y ~ output.u
    ]
    return ODESystem(ModelingToolkit.scalarize.(eqs), t, vcat(u..., y...), []; name, systems=[input, output])
end