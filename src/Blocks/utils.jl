@connector function RealInput(; name, nin = 1, u_start = nin > 1 ? zeros(nin) : 0.0)
    if nin == 1
        @variables u(t)=u_start [
            input = true,
            description = "Inner variable in RealInput $name"
        ]
    else
        @variables u(t)[1:nin]=u_start [
            input = true,
            description = "Inner variable in RealInput $name"
        ]
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name = name)
end
@doc """
    RealInput(;name, nin, u_start)

Connector with one input signal of type Real.

# Parameters:
- `nin=1`: Number of inputs
- `u_start=0`: Initial value for `u`

# States:
- `u`: Value of the connector; if nin=1 this is a scalar
""" RealInput

@connector function RealOutput(; name, nout = 1, u_start = nout > 1 ? zeros(nout) : 0.0)
    if nout == 1
        @variables u(t)=u_start [
            output = true,
            description = "Inner variable in RealOutput $name"
        ]
    else
        @variables u(t)[1:nout]=u_start [
            output = true,
            description = "Inner variable in RealOutput $name"
        ]
        u = collect(u)
    end
    ODESystem(Equation[], t, [u...], []; name = name)
end
@doc """
    RealOutput(;name, nout, u_start)

Connector with one output signal of type Real.

# Parameters:
- `nout=1`: Number of inputs
- `u_start=0`: Initial value for `u`

# States:
- `u`: Value of the connector; if nout=1 this is a scalar
""" RealOutput

"""
    SISO(;name, u_start = 0.0, y_start = 0.0)

Single input single output (SISO) continuous system block.

# Parameters:

  - `u`: Initial value for the input
  - `y`: Initial value for the output
"""
@mtkmodel SISO begin
    @parameters begin
        u_start = 0.0
        y_start = 0.0
    end
    @variables begin
        u(t) = u_start, [description = "Input of SISO system"]
        y(t) = y_start, [description = "Output of SISO system"]
    end
    @components begin
        input = RealInput(u_start = u_start)
        output = RealOutput(u_start = y_start)
    end
    @equations begin
        u ~ input.u
        y ~ output.u
    end
end

"""
    MIMO(; name, nin = 1, nout = 1, u_start = zeros(nin), y_start = zeros(nout))

Base class for a multiple input multiple output (MIMO) continuous system block.

# Parameters:

  - `nin`: Input dimension
  - `nout`: Output dimension
  - `u_start`: Initial value for the input
  - `y_start`: Initial value for the output
"""
@component function MIMO(; name, nin = 1, nout = 1, u_start = zeros(nin),
        y_start = zeros(nout))
    @named input = RealInput(nin = nin, u_start = u_start)
    @named output = RealOutput(nout = nout, u_start = y_start)
    @variables(u(t)[1:nin]=u_start, [description = "Input of MIMO system $name"],
        y(t)[1:nout]=y_start, [description = "Output of MIMO system $name"],)
    eqs = [
        [u[i] ~ input.u[i] for i in 1:nin]...,
        [y[i] ~ output.u[i] for i in 1:nout]...
    ]
    return ODESystem(eqs, t, vcat(u..., y...), []; name = name, systems = [input, output])
end
