@connector function RealInput(;name, nin=1, u_start=nin > 1 ? zeros(nin) : 0.0)
    if nin == 1
        @variables u(t) = u_start [input=true]
    else
        @variables u[1:nin](t) = u_start [input=true]
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name=name)
end
@doc """
    RealInput(;name, nin=1, u_start=nin > 1 ? 0.0 : zeros(nin))

Connector with one input signal of type Real.

# Parameters:
- `nin`: Number of inputs
- `u_start`: Initial value for `u`  

# States:
- `u`: Value of of the connector; if nin=1 this is a scalar
""" RealInput

@connector function RealOutput(;name, nout=1, u_start=nout > 1 ? zeros(nout) : 0.0)
    if nout == 1
        @variables u(t) = u_start [output=true]
    else
        @variables u[1:nout](t) = u_start [output=true]
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name=name)
end
@doc """
    RealOutput(;name, nout=1, u_start=nout > 1 ? 0.0 : zeros(nout))

Connector with one output signal of type Real.

# Parameters:
- `nout`: Number of inputs
- `u_start`: Initial value for `u`  

# States:
- `u`: Value of of the connector; if nout=1 this is a scalar
""" RealOutput

"""
    SISO(;name, u_start=0.0, y_start=0.0)

Single Input Single Output continuous control block.

# Parameters:
- `u_start`: Initial value for the input
- `y_start`: Initial value for the output
"""
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

"""
    MIMO(;name, nin=1, nout=1, u_start=zeros(nin), y_start=zeros(nout))

Base class for a multiple Input multiple Output continuous control block.

# Parameters:
- `nin`: Input dimension
- `nout`: Output dimension
- `u_start`: Initial value for the input
- `y_start`: Initial value for the output
"""
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
