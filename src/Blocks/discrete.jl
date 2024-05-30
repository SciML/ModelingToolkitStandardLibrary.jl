z = ShiftIndex()

"""
    DiscreteIntegrator(;name, k = 1, x = 0.0, method = :backward)

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
        method = :backward
    end
    @variables begin
        x(t) = 0.0, [description = "State of Integrator"]
    end
    @parameters begin
        k = 1, [description = "Gain"]
    end
    @structural_parameters begin
        Ts = ModelingToolkit.SampleTime()
        z = ShiftIndex()
    end
    @equations begin
        if method === :forward
            x(z) ~ x(z - 1) + k * Ts * u(z - 1)
        elseif method === :backward
            x(z) ~ x(z - 1) + k * Ts * u(z)
        elseif method ∈ (:trapezoidal, :tustin)
            x(z) ~ x(z - 1) + k * Ts * (u(z) + u(z - 1)) / 2
        end
        y ~ x(z)
    end
end

"""
    DiscreteDerivative(; k = 1)

A discrete-time derivative filter, corresponding to the transfer function
``k (z-1) / (T_s z)``

where `T_s` is the sample time of the derivative filter.

# Connectors:
- `input`
- `output`

# Parameters:
 - `k`: Gain of derivative filter
"""
@mtkmodel DiscreteDerivative begin
    @extend u, y = siso = SISO()
    @parameters begin
        k = 1, [description = "Gain"]
    end
    @structural_parameters begin
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @equations begin
        y(z) ~ k * (u(z) - u(z - 1)) / Ts
    end
end

"""
    Delay(; n = 1)

A discrete-time delay of `n` samples, corresponding to the transfer function ``z^{-n}``.

# Connectors:
- `input`
- `output`

# Structural Parameters:
 - `n`: Number of delay samples
"""
@mtkmodel Delay begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        n = 1
        z = ShiftIndex()
    end
    @equations begin
        y(z) ~ u(z - n)
    end
end

"""
    Difference()

A discrete-time difference operator, corresponding to the transfer function ``1 - z^{-1}``.

# Connectors:
- `input`
- `output`

For a discrete-time finite-difference derivative approximation, see [`DiscreteDerivative`](@ref).
"""
@mtkmodel Difference begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
    end
    @equations begin
        y(z) ~ u(z) - u(z - 1)
    end
end

"""
    ZeroOrderHold()

Zero-order-Hold translates a discrete time (clocked) signal into continuous time by holding the value of the discrete signal constant until the next sample.

# Connectors:
- `input` (discrete-time signal)
- `output` (continuous-time signal)
"""
@mtkmodel ZeroOrderHold begin
    @extend u, y = siso = SISO()
    @equations begin
        y ~ Hold(u)
    end
end

"""
    Sampler()
    Sampler(; dt::Real)
    Sampler(; clock::AbstractClock)

`Sampler` transforms a continuous-time signal into discrete time by sampling the input signal every time the associated clock ticks. The clock can be specified explicitly using the `clock` keyword argument, or implicitly by providing a sample time `dt`, in which case a standard periodic `Clock` is used. 

# Connectors:
- `input` (continuous-time signal)
- `output` (discrete-time signal)
"""
@mtkmodel Sampler begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        dt = nothing
        clock = (dt === nothing ? InferredDiscrete() : Clock(t, dt))
    end
    @equations begin
        y ~ Sample(clock)(u)
    end
end

function clock(s)
    ModelingToolkit.get_gui_metadata(s).type == GlobalRef(@__MODULE__, :Sampler) ||
        error("clock(sys) is supposed to be called on a ModelingToolkitStandardLibrary.Blocks.Sampler system, got $s")
    for eq in equations(s)
        td = ModelingToolkit.get_time_domain(eq.rhs)
        if td isa ModelingToolkit.AbstractClock
            return td
        end
    end
    error("No clock found")
end

"""
    ClockChanger()
    ClockChanger(; to::AbstractClock, from::AbstractClock)

# Connectors:
- `input` (continuous-time signal)
- `output` (discrete-time signal)
"""
@mtkmodel ClockChanger begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        to
        from
    end
    @equations begin
        y ~ ModelingToolkit.ClockChange(; to, from)(u)
    end
end

"""
    DiscretePIDParallel(;name, kp = 1, ki = 1, kd = 1, Ni = √(max(kd * ki, 1e-6)), Nd = 10kp, u_max = Inf, u_min = -u_max, wp = 1, wd = 1, Ts = 1, with_I = true, with_D = true, Imethod = :forward, Dmethod = :backward)

Discrete-time PID controller on parallel form with anti-windup and set-point weighting.

The controller is implemented on parallel form:

Simplified:
```math
u = k_p e + \\int k_i e dt + k_d \\dfrac{de}{dt} 
```

Detailed:
```math
u = k_p(w_p r - y) + \\int \\big( k_i (r - y) + N_i e_s \\big ) dt + k_d \\dfrac{d}{dt}(w_d r - y)
```

where `e_s = u - v` is the saturated error signal, `v` is the unsaturated control signal and `u` is the saturated control signal.

The derivative is filtered to allow a maximum gain of ``N_d``.

The integrator is discretized using the method specified by `Imethod`, options include
- `Imethod = :forward` (default): Corresponding to the transfer function ``T_s / (z - 1)``
- `Imethod = :backward`: Corresponding to the transfer function ``T_s z / (z - 1)``
- `Imethod = :trapezoidal`: Corresponding to the transfer function ``(T_s / 2) (z + 1) / (z - 1)``

The derivative is discretized using the method specified by `Dmethod`, options include
- `Dmethod = :forward`: Corresponding to the transfer function ``\\dfrac{N (z-1)}{z - \\dfrac{k_d-N T_s}{k_d}}``.
- `Dmethod = :backward` (default): Corresponding to the transfer function ``\\dfrac{\\dfrac{Nk_d}{k_d + N T_s}(z-1)}{z - \\dfrac{k_d}{k_d + N T_s}}``

Anti windup is realized by tracking using the gain ``N_i`` on the error signal ``e_s`` when the output is saturated.

To use the controller in 1DOF mode, i.e., with only the control error as input, connect the error signal to the `reference` connector, connect a `Constant(; k = 0)` to the `measurement` connector and set `wp = wd = 1`.

# Connectors:
- `reference`: The reference signal to the controller (or the error signal if used in 1DOF mode)
- `measurement`: The measurement feedback
- `ctr_output`: The control signal output

# Parameters:
- `kp`: Proportional gain
- `ki`: Integral gain (only active if `with_I = true`)
- `kd`: Derivative gain (only active if `with_D = true`)
- `Ni`: Anti-windup gain (only active if `with_I = true`)
- `Nd`: Maximum derivative gain (only active if `with_D = true`). Typically set to 10-100 times the proportional gain.
- `u_max`: Maximum output above which the output is saturated
- `u_min`: Minimum output below which the output is saturated. This defaults to `-u_max` if `u_max > 0` and `-Inf` otherwise.
- `wp`: `[0, 1]` Set-point weighting in the proportional part. Set to `0` to prevent step changes in the output due to step changes in the reference.
- `wd`: `[0, 1]` Set-point weighting in the derivative part. Set to `0` to prevent very large impulsive changes in the output due to step changes in the reference.
- `with_I`: Whether or not to include the integral part
- `with_D`: Whether or not to include the derivative part
- `Imethod`: Discretization method for the integrator (see details above)
- `Dmethod`: Discretization method for the derivative (see details above)


# Extended help:
## Internal variables:
- `I`: State of integrator
- `D`: State of filtered derivative
- `r`: Reference signal internal variable
- `y`: Measurement signal internal variable
- `wde`: Setpoint-weighted error for derivative
- `v`: Un-saturated output of the controller
- `u`: Saturated output of the controller
- `eI`: Error signal input to integrator including anit-windup tracking signal
- `e`: Error signal
"""
@mtkmodel DiscretePIDParallel begin
    @structural_parameters begin
        Imethod = :forward
        Dmethod = :backward
        with_I = true
        with_D = true
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @components begin
        reference = RealInput()
        measurement = RealInput()
        ctr_output = RealOutput()
    end
    @variables begin
        I(t) = 0.0, [description = "State of Integrator"]
        D(t) = 0.0, [description = "State of filtered derivative"]
        r(t), [guess = 0, description = "Reference signal internal variable"]
        y(t), [guess = 0, description = "Measurement signal internal variable"]
        wde(t), [guess = 0, description = "Setpoint-weighted error for derivative"]
        v(t), [guess = 0, description = "Un-saturated output of the controller"]
        u(t), [guess = 0, description = "Saturated output of the controller"]
        eI(t),
        [guess = 0,
            description = "Error signal input to integrator including anit-windup tracking signal"]
        e(t), [guess = 0, description = "Error signal"]
    end
    @parameters begin
        kp = 1, [description = "Proportional gain"]
        ki = 1, [description = "Integral gain"]
        kd = 1, [description = "Derivative gain"]
        Ni = √(max(kd * ki, 1e-6)), [description = "Anti-windup gain"]
        Nd = 10 * kp, [description = "Maximum derivative gain"]
        u_max = Inf, [description = "Maximum output"]
        u_min = ifelse(u_max > 0, -u_max, -Inf), [description = "Minimum output"]
        wp = 1, [description = "Set-point weighting in the proportional part."]
        wd = 1, [description = "Set-point weighting in the derivative part."]
    end
    @equations begin
        r ~ reference.u
        y ~ measurement.u
        u ~ ctr_output.u
        e ~ r - y
        v ~ kp * (wp * r - y) + I(z - 1) + D # Unsaturated control signal
        u ~ _clamp(v, u_min, u_max) # Saturated control signal
        if with_I
            eI ~ e + Ni * (u - v) # Add anti-windup tracking signal to error before integration
            if Imethod === :forward
                I(z) ~ I(z - 1) + Ts * ki * eI(z - 1)
            elseif Imethod === :backward
                I(z) ~ I(z - 1) + Ts * ki * eI(z)
            elseif Imethod ∈ (:trapezoidal, :tustin)
                I(z) ~ I(z - 1) + Ts * ki * (eI(z) + eI(z - 1)) / 2
            else
                error("Unknown integrator discretization method $Imethod, must be one of :forward, :backward, :trapezoidal")
            end
        else
            I(z) ~ 0
        end
        if with_D
            wde ~ wd * r - y
            if Dmethod === :forward
                D(z) ~ (kd - Nd * Ts) / kd * D(z - 1) + Nd * (wde(z) - wde(z - 1))
            elseif Dmethod === :backward
                D(z) ~ kd / (kd + Nd * Ts) * D(z - 1) +
                       Nd * kd / (kd + Nd * Ts) * (wde(z) - wde(z - 1))
            else
                error("Unknown derivative discretization method $Dmethod, must be one of :forward, :backward")
            end
        else
            D(z) ~ 0
        end
    end
end

"""
    DiscretePIDStandard(;name, K = 1, Ti = 1, Td = 1, Ni = √(max(kd * ki, 1e-6)), Nd = 10, u_max = Inf, u_min = -u_max, wp = 1, wd = 1, Ts = 1, with_I = true, with_D = true, Imethod = :forward, Dmethod = :backward)

Discrete-time PID controller on standard (ideal) form with anti-windup and set-point weighting.

The controller is implemented on standard form

Simplified:
```math
u = K \\left( e + \\int \\frac{1}{T_i} e  dt + T_d \\dfrac{de}{dt} \\right)
```

Detailed:
```math
u = K \\left( (w_p r - y) + \\int \\big( \\frac{1}{T_i} (r - y) + N_i e_s \\big ) dt + T_d \\dfrac{d}{dt}(w_d r - y) \\right)
```

where `e_s = u - v` is the saturated error signal, `v` is the unsaturated control signal and `u` is the saturated control signal.

The derivative is filtered to allow a maximum gain of ``N_d``.

The integrator is discretized using the method specified by `Imethod`, options include
- `Imethod = :forward` (default): Corresponding to the transfer function ``T_s / (z - 1)``
- `Imethod = :backward`: Corresponding to the transfer function ``T_s z / (z - 1)``
- `Imethod = :trapezoidal`: Corresponding to the transfer function ``(T_s / 2) (z + 1) / (z - 1)``

The derivative is discretized using the method specified by `Dmethod`, options include
- `Dmethod = :forward`: Corresponding to the transfer function ``\\dfrac{N (z-1)}{z - \\dfrac{T_d-N T_s}{T_d}}``.
- `Dmethod = :backward` (default): Corresponding to the transfer function ``\\dfrac{\\dfrac{NT_d}{T_d + N T_s}(z-1)}{z - \\dfrac{T_d}{T_d + N T_s}}``

Anti windup is realized by tracking using the gain ``N_i`` on the error signal ``e_s`` when the output is saturated.

To use the controller in 1DOF mode, i.e., with only the control error as input, connect the error signal to the `reference` connector, connect a `Constant(; k = 0)` to the `measurement` connector and set `wp = wd = 1`.

# Connectors:
- `reference`: The reference signal to the controller (or the error signal if used in 1DOF mode)
- `measurement`: The measurement feedback
- `ctr_output`: The control signal output

# Parameters:
- `K`: Proportional gain
- `Ti`: Integral time constant (only active if `with_I = true`)
- `Td`: Derivative time (only active if `with_D = true`)
- `Ni`: Anti-windup gain (only active if `with_I = true`)
- `Nd`: Maximum derivative gain (only active if `with_D = true`). Typically set to 10-100.
- `u_max`: Maximum output above which the output is saturated
- `u_min`: Minimum output below which the output is saturated. This defaults to `-u_max` if `u_max > 0` and `-Inf` otherwise.
- `wp`: `[0, 1]` Set-point weighting in the proportional part. Set to `0` to prevent step changes in the output due to step changes in the reference.
- `wd`: `[0, 1]` Set-point weighting in the derivative part. Set to `0` to prevent very large impulsive changes in the output due to step changes in the reference.
- `with_I`: Whether or not to include the integral part
- `with_D`: Whether or not to include the derivative part
- `Imethod`: Discretization method for the integrator (see details above)
- `Dmethod`: Discretization method for the derivative (see details above)


# Extended help:
## Internal variables:
- `I`: State of integrator
- `D`: State of filtered derivative
- `r`: Reference signal internal variable
- `y`: Measurement signal internal variable
- `wde`: Setpoint-weighted error for derivative
- `v`: Un-saturated output of the controller
- `u`: Saturated output of the controller
- `eI`: Error signal input to integrator including anit-windup tracking signal
- `e`: Error signal
"""
@mtkmodel DiscretePIDStandard begin
    @structural_parameters begin
        Imethod = :forward
        Dmethod = :backward
        with_I = true
        with_D = true
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @components begin
        reference = RealInput()
        measurement = RealInput()
        ctr_output = RealOutput()
    end
    @variables begin
        I(t) = 0.0, [description = "State of Integrator"]
        D(t) = 0.0, [description = "State of filtered derivative"]
        r(t) = 0.0, [description = "Reference signal internal variable"]
        y(t) = 0.0, [description = "Measurement signal internal variable"]
        wde(t) = 0.0, [description = "Setpoint-weighted error for derivative"]
        v(t) = 0.0, [description = "Un-saturated output of the controller"]
        u(t) = 0.0, [description = "Saturated output of the controller"]
        eI(t) = 0.0,
        [description = "Error signal input to integrator including anit-windup tracking signal"]
        e(t) = 0.0, [description = "Error signal"]
    end
    @parameters begin
        K = 1, [description = "Proportional gain"]
        Ti = 1, [description = "Integral time"]
        Td = 1, [description = "Derivative time"]
        Ni = √(max(Td / Ti, 1e-6)), [description = "Anti-windup gain"]
        Nd = 10, [description = "Maximum derivative gain"]
        u_max = Inf, [description = "Maximum output"]
        u_min = ifelse(u_max > 0, -u_max, -Inf), [description = "Minimum output"]
        wp = 1, [description = "Set-point weighting in the proportional part."]
        wd = 1, [description = "Set-point weighting in the derivative part."]
    end
    @equations begin
        r ~ reference.u
        y ~ measurement.u
        u ~ ctr_output.u
        e ~ r - y
        v ~ K * ((wp * r - y) + I(z - 1) + D) # Unsaturated control signal
        u ~ _clamp(v, u_min, u_max) # Saturated control signal
        if with_I
            eI ~ e + Ni * (u - v) # Add anti-windup tracking signal to error before integration
            if Imethod === :forward
                I(z) ~ I(z - 1) + Ts / Ti * eI(z - 1)
            elseif Imethod === :backward
                I(z) ~ I(z - 1) + Ts / Ti * eI(z)
            elseif Imethod ∈ (:trapezoidal, :tustin)
                I(z) ~ I(z - 1) + Ts / Ti * (eI(z) + eI(z - 1)) / 2
            else
                error("Unknown integrator discretization method $Imethod, must be one of :forward, :backward, :trapezoidal")
            end
        else
            I(z) ~ 0
        end
        if with_D
            wde = wd * r - y
            if Dmethod === :forward
                D(z) ~ (Td - Nd * Ts) / Td * D(z - 1) + Nd * (wde(z) - wde(z - 1))
            elseif Dmethod === :backward
                D(z) ~ Td / (Td + Nd * Ts) * D(z - 1) +
                       Nd * Td / (Td + Nd * Ts) * (wde(z) - wde(z - 1))
            else
                error("Unknown derivative discretization method $Dmethod, must be one of :forward, :backward")
            end
        else
            D(z) ~ 0
        end
    end
end

# """
#     DiscreteStateSpace(; A, B, C, D, u0 = zeros(size(B, 2)), y0 = zeros(size(C, 1)))

# A linear, time invariant, discrete-time state-space system on the form
# ```math
# x_{k+1} = A x_k + B u_k
# y_k = C x_k + D u_k
# ```

# `y0` and `u0` can be used to set an operating point, providing them changes the dynamics from an LTI system to the affine system
# ```math
# x_{k+1} = A x_k + B (u_k - u_0)
# y_k = C x_k + D (u_k - u_0) + y_0
# ```
# """
# @mtkmodel DiscreteStateSpace begin
#     @structural_parameters begin
#         z = ShiftIndex()
#     end
#     @parameters begin
#         A, [description = "State matrix"]
#         B, [description = "Input matrix"]
#         C, [description = "Output matrix"]
#         (D = 0),      [description = "Feedthrough matrix"]
#         (u0[1:size(B, 2)] = 0), [description = "Input operating point"]
#         (y0[1:size(C, 1)] = 0), [description = "Output operating point"]
#     end
#     begin
#         nx = size(A, 1)
#         nu = size(B, 2)
#         ny = size(C, 1)
#         size(A, 2) == nx || error("`A` has to be a square matrix.")
#         size(B, 1) == nx || error("`B` has to be of dimension ($nx x $nu).")
#         size(C, 2) == nx || error("`C` has to be of dimension ($ny x $nx).")
#     end
#     @components begin
#         input = RealInput(nin = nu)
#         output = RealOutput(nout = ny)
#     end
#     @variables begin
#         (x(t)[1:nx]), [description = "State"]
#         (u(t)[1:nu]), [description = "Input"]
#         (y(t)[1:ny]), [description = "Output"]
#     end
#     @equations begin
#         x(z) ~ A * x(z-1) .+ B * (u(z-1) .- u0)
#         y(z) ~ C * x(z) .+ D * (u(z) .- u0) .+ y0
#         output.u ~ y
#         input.u ~ u
#     end
# end

# function DiscreteStateSpace(A, B, C, D = nothing; kwargs...)
#     DiscreteStateSpace(; A, B, C, D, kwargs...)
# end

"""
    DiscreteTransferFunction(; b, a)

A discrete-time transfer function on the form
```math
H(z) = \\dfrac{B(z)}{A(z)} = \\dfrac{b_{n_b} z^{n_b} + b_{n_b-1} z^{n_b-1} + \\cdots + b_1 z + b_0}{a_{n_a} z^{n_a} + a_{n_a-1} z^{n_a-1} + \\cdots + a_1 z + a_0}
```

With the coeffiencents specified in decreasing orders of ``z``, i.e., ``b = [b_{n_b}, b_{n_b-1}, \\cdots, b_1, b_0]`` and ``a = [a_{n_a}, a_{n_a-1}, \\cdots, a_1, a_0]``.

## Parameters:
- `b`: Numerator coefficients of transfer function (e.g., z - 1 is specified as `[1,-1]`)
- `a`: Denominator coefficients of transfer function (e.g., z^2 - 0.78z + 0.37 is specified as `[1, -0.78, 0.37]`)
- `verbose`: Whether or not to print a warning if the numerator degree is larger than the denominator degree.

## Connectors:
- `input`: Input signal
- `output`: Output signal

# Extended help:
This component supports SISO systems only. To simulate MIMO transfer functions, use [ControlSystemsBase.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/man/creating_systems/) to convert the transfer function to a statespace system, optionally compute a minimal realization using [`minreal`](https://juliacontrol.github.io/ControlSystems.jl/stable/lib/constructors/#ControlSystemsBase.minreal), and then use a [`DiscreteStateSpace`](@ref) component instead.

See also [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/stable/) for an interface between [ControlSystems.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/) and ModelingToolkit.jl for advanced manipulation of transfer functions and linear statespace systems. For linearization, see [`linearize`](@ref) and [Linear Analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/).
"""
# @mtkmodel DiscreteTransferFunction begin
#     @parameters begin
#         b = [1], [description = "Numerator coefficients of transfer function (e.g., z - 1 is specified as [1,-1])"]
#         a = [1], [description = "Denominator coefficients of transfer function (e.g., z^2 - 0.78z + 0.37 is specified as [1, -0.78, 0.37])"]
#     end
#     @structural_parameters begin
#         verbose = true
# Ts = SampleTime()
# z = ShiftIndex()
#     end
#     begin
#         na = length(a)
#         nb = length(b)
#         verbose && nb > na && @warn "DiscreteTransferFunction: Numerator degree is larger than denominator degree, this is not a proper transfer function. Simulation of a model including this transfer funciton will require at least $(nb-na) samples additional delay in series. Silence this warning with verbose=false"
#         Ts = SampleTime()
#     end
#     @components begin
#         input = RealInput()
#         output = RealOutput()
#     end
#     @variables begin
#         u(t) = 0.0, [description = "Input flowing through connector `input`"]
#         y(t) = 0.0, [description = "Output flowing through connector `output`"]
#     end
#     @equations begin
#         sum(y(z+k-1) * a[na-k+1] for k in 1:na) ~ sum(u(z+k-1) * b[nb-k+1] for k in 1:nb)
#         input.u ~ u
#         output.u ~ y
#     end
# end

# DiscreteTransferFunction(b, a; kwargs...) = DiscreteTransferFunction(; b, a, kwargs...)

##
