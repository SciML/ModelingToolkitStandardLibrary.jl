@connector function RealInput(;name, nin=1, u_start=nin > 1 ? 0.0 : zeros(nin))
    if nin == 1
        @variables u(t) = u_start
    else
        @variables u[1:nin](t) = u_start
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name=name)
end

@connector function RealOutput(;name, nout=1, u_start=nout > 1 ? 0.0 : zeros(nout))
    if nout == 1
        @variables u(t) = u_start
    else
        @variables u[1:nout](t) = u_start
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name=name)
end

function SISO(;name, u_start=0.0, y_start=0.0)
    @named input = RealInput(u_start=u_start)
    @named output = RealOutput(u_start=y_start)
    @variables begin
        u(t)=u_start
        y(t)=y_start
    end
    eqs = [
        u ~ input.u
        y ~ output.u
    ]
    return ODESystem(eqs, t, [u, y], []; name=name, systems=[input, output])
end

function MIMO(;name, nin=1, nout=1, u_start=zeros(nin), y_start=zeros(nout))
    @named input = RealInput(nin=nin, u_start=u_start)
    @named output = RealOutput(nout=nout, u_start=y_start)
    @variables begin
        u[1:nin](t)=u_start
        y[1:nout](t)=y_start
    end
    eqs = [
        [u[i] ~ input.u[i] for i in 1:nin]...,
        [y[i] ~ output.u[i] for i in 1:nout]...,
    ]
    return ODESystem(eqs, t, vcat(u..., y...), []; name=name, systems=[input, output])
end