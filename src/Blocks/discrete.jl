z = ShiftIndex()

"""
    Integrator(;name, k = 1, x = 0.0, method = :forward)

Outputs `y = ∫k*u dt`, corresponding to the discrete-time transfer function
- `method = :forward`: ``T_s / (z - 1)``
- `method = :backward`: ``T_s z / (z - 1)``
- `method = :trapezoidal`: ``(T_s / 2) (z + 1) / (z - 1)``

where `T_s` is the sample time of the integrator.

Initial value of integrator state ``x`` can be set with `x`

# Connectors:

  - `input`
  - `output`

# Parameters:

  - `k`: Gain of integrator
"""
@mtkmodel DiscreteIntegrator begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        method = :forward
    end
    @variables begin
        x(t) = 0.0, [description = "State of Integrator"]
    end
    @parameters begin
        k = 1, [description = "Gain"]
    end
    begin
        Ts = sampletime(x)
    end
    @equations begin
        if method === :forward
            x(z) ~ x(z-1) + k * Ts * u(z-1)
        elseif method === :backward
            x(z) ~ x(z-1) + k * Ts * u(z)
        elseif method === :trapezoidal
            x(z) ~ x(z-1) + k * Ts * (u(z) + u(z-1)) / 2
        end
        y ~ x(z)
    end
end

@mtkmodel DiscreteDerivative begin
    @extend u, y = siso = SISO()
    @parameters begin
        k = 1, [description = "Gain"]
    end
    begin
        Ts = sampletime(u)
    end
    @equations begin
        y(z) ~ k*(u(z) - u(z-1)) / Ts
    end
end


@mtkmodel Delay begin
    @extend u, y = siso = SISO()
    @parameters begin
        n = 1, [description = "Number of delay samples"]
    end
    @equations begin
        y ~ u(z-n)
    end
end

@mtkmodel Difference begin
    @extend u, y = siso = SISO()
    @equations begin
        y(z) ~ u(z) - u(z-1)
    end
end


@mtkmodel ZeroOrderHold begin
    @extend u, y = siso = SISO()
    @equations begin
        y ~ Hold(u)
    end
end

@mtkmodel Sampler begin
    @extend u, y = siso = SISO()
    # @parameters begin
    #     Ts = 1, [description = "Sample interval"]
    # end # TODO: figure out how to connect a clock
    @equations begin
        y ~ Sample(t)(u)
    end
end
