"""
    Integrator(;name, k=1, x_start=0.0)

Outputs `y = ∫k*u dt`, corresponding to the transfer function `1/s`.

# Connectors:
- `input`
- `output`

# Parameters:
- `k`: Gain of integrator
- `x_start`: Initial value of integrator
"""
function Integrator(; name, k = 1, x_start = 0.0)
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t) = x_start
    pars = @parameters k = k
    eqs = [D(x) ~ k * u
           y ~ x]
    extend(ODESystem(eqs, t, sts, pars; name = name), siso)
end

"""
    Derivative(; name, k=1, T, x_start=0.0)

Outputs an approximate derivative of the input. The transfer function of this block is
```
k       k     
─ - ──────────
T    2 ⎛    1⎞
    T ⋅⎜s + ─⎟
       ⎝    T⎠
```
and a state-space realization is given by `ss(-1/T, 1/T, -k/T, k/T)`
where `T` is the time constant of the filter.
A smaller `T` leads to a more ideal approximation of the derivative.

# Parameters:
- `k`: Gain
- `T`: [s] Time constants (T>0 required; T=0 is ideal derivative block)
- `x_start`: Initial value of state

# Connectors:
- `input`
- `output`
"""
function Derivative(; name, k = 1, T, x_start = 0.0)
    T > 0 || throw(ArgumentError("Time constant `T` has to be strictly positive"))
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t) = x_start
    pars = @parameters T=T k=k
    eqs = [D(x) ~ (u - x) / T
           y ~ (k / T) * (u - x)]
    extend(ODESystem(eqs, t, sts, pars; name = name), siso)
end

"""
    FirstOrder(; name, k=1, T, x_start=0.0, lowpass=true)

A first-order filter with a single real pole in `s = -T` and gain `k`. If `lowpass=true` (default), the transfer function
is given by `Y(s)/U(s) = `
```
   k   
───────
sT + 1
```
and if `lowpass=false`, by
```
sT + 1 - k   
──────────
  sT + 1
```


# Parameters:
- `k`: Gain
- `T`: [s] Time constants (T>0 required)
- `x_start`: Initial value of state

# Connectors:
- `input`
- `output`

See also [`SecondOrder`](@ref)
"""
function FirstOrder(; name, k = 1, T, x_start = 0.0, lowpass = true)
    T > 0 || throw(ArgumentError("Time constant `T` has to be strictly positive"))
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t) = x_start
    pars = @parameters T=T k=k
    eqs = [D(x) ~ (k * u - x) / T
           lowpass ? y ~ x : y ~ k * u - x]
    extend(ODESystem(eqs, t, sts, pars; name = name), siso)
end

"""
    SecondOrder(; name, k=1, w, d, x_start=0.0, xd_start=0.0)

A second-order filter with gain `k`, a bandwidth of `w` rad/s and relative damping `d`. The transfer function
is given by `Y(s)/U(s) = `
```
      k*w^2   
─────────────────
s² + 2d*w*s + w^2
```
Critical damping corresponds to `d=1`, which yields the fastest step response without overshoot, d < 1` results in an under-damped filter while `d > 1` results in an over-damped filter.
`d = 1/√2` corresponds to a Butterworth filter of order 2 (maximally flat frequency response).

# Parameters:
- `k`: Gain
- `w`: [`rad/s`] Angular frequency
- `d`: Damping
- `x_start`: Initial value of state (output)
- `xd_start`: Initial value of derivative of state (output)

# Connectors:
- `input`
- `output`
"""
function SecondOrder(; name, k = 1, w, d, x_start = 0.0, xd_start = 0.0)
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t)=x_start xd(t)=xd_start
    pars = @parameters k=k w=w d=d
    eqs = [D(x) ~ xd
           D(xd) ~ w * (w * (k * u - x) - 2 * d * xd)
           y ~ x]
    extend(ODESystem(eqs, t, sts, pars; name = name), siso)
end

"""
    PI(;name, k=1, T, x_start=0.0)

Textbook version of a PI-controller without actuator saturation and anti-windup measure.

# Parameters:
- `k`: Gain
- `T`: [s] Integrator time constant (T>0 required)
- `x_start`: Initial value for the integrator

# Connectors:
- `err_input`
- `ctr_output`

See also [`LimPI`](@ref)
"""
function PI(; name, k = 1, T, x_start = 0.0)
    T > 0 || throw(ArgumentError("Time constant `T` has to be strictly positive"))
    @named err_input = RealInput() # control error
    @named ctr_output = RealOutput() # control signal
    @named gainPI = Gain(k)
    @named addPI = Add()
    @named int = Integrator(k = 1 / T, x_start = x_start)
    sys = [err_input, ctr_output, gainPI, addPI, int]
    eqs = [
        connect(err_input, addPI.input1),
        connect(addPI.output, gainPI.input),
        connect(gainPI.output, ctr_output),
        connect(err_input, int.input),
        connect(int.output, addPI.input2),
    ]
    ODESystem(eqs, t, [], []; name = name, systems = sys)
end

"""
    PID(;name, k=1, Ti=false, Td=false, Nd=10, xi_start=0, xd_start=0)

Text-book version of a PID-controller without actuator saturation and anti-windup measure.

# Parameters:
- `k`: Gain
- `Ti`: [s] Integrator time constant (Ti>0 required). If set to false no integral action is used.
- `Td`: [s] Derivative time constant (Td>0 required). If set to false no derivative action is used.
- `Nd`: [s] Time constant for the derivative approximation (Nd>0 required; Nd=0 is ideal derivative).
- `x_start`: Initial value for the integrator.
- `xd_start`: Initial value for the derivative state.

# Connectors:
- `err_input`
- `ctr_output`

See also [`LimPID`](@ref)
"""
function PID(; name, k = 1, Ti = false, Td = false, Nd = 10, xi_start = 0, xd_start = 0)
    with_I = !isequal(Ti, false)
    with_D = !isequal(Td, false)
    @named err_input = RealInput() # control error
    @named ctr_output = RealOutput() # control signal
    !isequal(Ti, false) &&
        (Ti ≥ 0 || throw(ArgumentError("Ti out of bounds, got $(Ti) but expected Ti ≥ 0")))
    !isequal(Td, false) &&
        (Td ≥ 0 || throw(ArgumentError("Td out of bounds, got $(Td) but expected Td ≥ 0")))
    Nd > 0 || throw(ArgumentError("Nd out of bounds, got $(Nd) but expected Nd > 0"))

    @named gainPID = Gain(k)
    @named addPID = Add3()
    if with_I
        @named int = Integrator(k = 1 / Ti, x_start = xi_start)
    else
        @named Izero = Constant(k = 0)
    end
    if with_D
        @named der = Derivative(k = 1 / Td, T = 1 / Nd, x_start = xd_start)
    else
        @named Dzero = Constant(k = 0)
    end
    sys = [err_input, ctr_output, gainPID, addPID]
    if with_I
        push!(sys, int)
    else
        push!(sys, Izero)
    end
    if with_D
        push!(sys, der)
    else
        push!(sys, Dzero)
    end
    eqs = [
        connect(err_input, addPID.input1),
        connect(addPID.output, gainPID.input),
        connect(gainPID.output, ctr_output),
    ]
    if with_I
        push!(eqs, connect(err_input, int.input))
        push!(eqs, connect(int.output, addPID.input2))
    else
        push!(eqs, connect(Izero.output, addPID.input2))
    end
    if with_D
        push!(eqs, connect(err_input, der.input))
        push!(eqs, connect(der.output, addPID.input3))
    else
        push!(eqs, connect(Dzero.output, addPID.input3))
    end
    ODESystem(eqs, t, [], []; name = name, systems = sys)
end

"""
    LimPI(;name, k=1, T, u_max=1, u_min=-u_max, Ta)

Text-book version of a PI-controller with actuator saturation and anti-windup measure.

# Parameters:
- `k`: Gain
- `T`: [s] Integrator time constant (T>0 required)
- `Ta`: [s] Tracking time constant (Ta>0 required)
- `x_start`: Initial value for the integrator

# Connectors:
- `err_input`
- `ctr_output`
"""
function LimPI(; name, k = 1, T, u_max, u_min = -u_max, Ta, x_start = 0.0)
    Ta > 0 || throw(ArgumentError("Time constant `Ta` has to be strictly positive"))
    T > 0 || throw(ArgumentError("Time constant `T` has to be strictly positive"))
    u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))
    @named err_input = RealInput() # control error
    @named ctr_output = RealOutput() # control signal
    @named gainPI = Gain(k)
    @named addPI = Add()
    @named addTrack = Add()
    @named int = Integrator(k = 1 / T, x_start = x_start)
    @named limiter = Limiter(y_max = u_max, y_min = u_min)
    @named addSat = Add(k1 = 1, k2 = -1)
    @named gainTrack = Gain(1 / Ta)
    sys = [err_input, ctr_output, gainPI, addPI, int, addTrack, limiter, addSat, gainTrack]
    eqs = [
        connect(err_input, addPI.input1),
        connect(addPI.output, gainPI.input),
        connect(gainPI.output, limiter.input),
        connect(limiter.output, ctr_output),
        connect(limiter.input, addSat.input2),
        connect(limiter.output, addSat.input1),
        connect(addSat.output, gainTrack.input),
        connect(err_input, addTrack.input1),
        connect(gainTrack.output, addTrack.input2),
        connect(addTrack.output, int.input),
        connect(int.output, addPI.input2),
    ]
    ODESystem(eqs, t, [], []; name = name, systems = sys)
end

"""
    LimPID(; k, Ti=false, Td=false, wp=1, wd=1, Ni, Nd=12, u_max=Inf, u_min=-u_max, gains = false, name)

Proportional-Integral-Derivative (PID) controller with output saturation, set-point weighting and integrator anti-windup.

The equation for the control signal is roughly
```
k(ep + 1/Ti * ∫e + 1/Td * d/dt(ed))
e = u_r - u_y
ep = wp*u_r - u_y
ed = wd*u_r - u_y
```
where the transfer function for the derivative includes additional filtering, see `? Derivative` for more details.

# Parameters:
- `k`: Proportional gain
- `Ti`: [s] Integrator time constant. Set to `false` to turn off integral action.
- `Td`: [s] Derivative time constant. Set to `false` to turn off derivative action.
- `wp`: [0,1] Set-point weighting in the proportional part.
- `wd`: [0,1] Set-point weighting in the derivative part.
- `Nd`: [1/s] Derivative limit, limits the derivative gain to Nd/Td. Reasonable values are ∈ [8, 20]. A higher value gives a better approximation of an ideal derivative at the expense of higher noise amplification.
- `Ni`: `Ni*Ti` controls the time constant `Ta` of anti-windup tracking. A common (default) choice is `Ta = √(Ti*Td)` which is realized by `Ni = √(Td / Ti)`. Anti-windup can be effectively turned off by setting `Ni = Inf`.
- `gains`: If `gains = true`, `Ti` and `Td` will be interpreted as gains with a fundamental PID transfer function on parallel form `ki=Ti, kd=Td, k + ki/s + kd*s`.

# Connectors:
- `reference`
- `measurement`
- `ctr_output`
"""
function LimPID(; name, k = 1, Ti = false, Td = false, wp = 1, wd = 1,
                Ni = Ti == 0 ? Inf : √(max(Td / Ti, 1e-6)),
                Nd = 10,
                u_max = Inf,
                u_min = u_max > 0 ? -u_max : -Inf,
                gains = false,
                xi_start = 0.0,
                xd_start = 0.0)
    with_I = !isequal(Ti, false)
    with_D = !isequal(Td, false)
    with_AWM = Ni != Inf
    if gains
        Ti = k / Ti
        Td = Td / k
    end
    0 ≤ wp ≤ 1 ||
        throw(ArgumentError("wp out of bounds, got $(wp) but expected wp ∈ [0, 1]"))
    0 ≤ wd ≤ 1 ||
        throw(ArgumentError("wd out of bounds, got $(wd) but expected wd ∈ [0, 1]"))
    !isequal(Ti, false) &&
        (Ti ≥ 0 || throw(ArgumentError("Ti out of bounds, got $(Ti) but expected Ti ≥ 0")))
    !isequal(Td, false) &&
        (Td ≥ 0 || throw(ArgumentError("Td out of bounds, got $(Td) but expected Td ≥ 0")))
    u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))
    Nd > 0 || throw(ArgumentError("Nd out of bounds, got $(Nd) but expected Nd > 0"))

    @named reference = RealInput()
    @named measurement = RealInput()
    @named ctr_output = RealOutput() # control signal
    @named addP = Add(k1 = wp, k2 = -1)
    @named gainPID = Gain(k)
    @named addPID = Add3()
    @named limiter = Limiter(y_max = u_max, y_min = u_min)
    if with_I
        if with_AWM
            @named addI = Add3(k1 = 1, k2 = -1, k3 = 1)
            @named addSat = Add(k1 = 1, k2 = -1)
            @named gainTrack = Gain(1 / (k * Ni))
        else
            @named addI = Add(k1 = 1, k2 = -1)
        end
        @named int = Integrator(k = 1 / Ti, x_start = xi_start)
    else
        @named Izero = Constant(k = 0)
    end
    if with_D
        @named der = Derivative(k = 1 / Td, T = 1 / Nd, x_start = xd_start)
        @named addD = Add(k1 = wd, k2 = -1)
    else
        @named Dzero = Constant(k = 0)
    end

    sys = [reference, measurement, ctr_output, addP, gainPID, addPID, limiter]
    if with_I
        if with_AWM
            push!(sys, [addSat, gainTrack]...)
        end
        push!(sys, [addI, int]...)
    else
        push!(sys, Izero)
    end
    if with_D
        push!(sys, [addD, der]...)
    else
        push!(sys, Dzero)
    end

    eqs = [
        connect(reference, addP.input1),
        connect(measurement, addP.input2),
        connect(addP.output, addPID.input1),
        connect(addPID.output, gainPID.input),
        connect(gainPID.output, limiter.input),
        connect(limiter.output, ctr_output),
    ]
    if with_I
        push!(eqs, connect(reference, addI.input1))
        push!(eqs, connect(measurement, addI.input2))
        if with_AWM
            push!(eqs, connect(limiter.input, addSat.input2))
            push!(eqs, connect(limiter.output, addSat.input1))
            push!(eqs, connect(addSat.output, gainTrack.input))
            push!(eqs, connect(gainTrack.output, addI.input3))
        end
        push!(eqs, connect(addI.output, int.input))
        push!(eqs, connect(int.output, addPID.input3))
    else
        push!(eqs, connect(Izero.output, addPID.input3))
    end
    if with_D
        push!(eqs, connect(reference, addD.input1))
        push!(eqs, connect(measurement, addD.input2))
        push!(eqs, connect(addD.output, der.input))
        push!(eqs, connect(der.output, addPID.input2))
    else
        push!(eqs, connect(Dzero.output, addPID.input2))
    end

    ODESystem(eqs, t, [], []; name = name, systems = sys)
end

"""
    StateSpace(A, B, C, D=0; x_start=zeros(size(A,1)), name)

A linear, time-invariant state-space system on the form.
```
ẋ = Ax + Bu
y = Cx + Du
```
Transfer functions can also be simulated by converting them to a StateSpace form.
"""
function StateSpace(; A, B, C, D = nothing, x_start = zeros(size(A, 1)), name)
    nx, nu, ny = size(A, 1), size(B, 2), size(C, 1)
    size(A, 2) == nx || error("`A` has to be a square matrix.")
    size(B, 1) == nx || error("`B` has to be of dimension ($nx x $nu).")
    size(C, 2) == nx || error("`C` has to be of dimension ($ny x $nx).")
    if B isa AbstractVector
        B = reshape(B, length(B), 1)
    end
    if isnothing(D) || iszero(D)
        D = zeros(ny, nu)
    else
        size(D) == (ny, nu) || error("`D` has to be of dimension ($ny x $nu).")
    end
    @named input = RealInput(nin = nu)
    @named output = RealOutput(nout = ny)
    @variables x(t)[1:nx] = x_start
    # pars = @parameters A=A B=B C=C D=D # This is buggy
    eqs = [ # FIXME: if array equations work
        [Differential(t)(x[i]) ~ sum(A[i, k] * x[k] for k in 1:nx) +
                                 sum(B[i, j] * input.u[j] for j in 1:nu) for i in 1:nx]..., # cannot use D here
        [output.u[j] ~ sum(C[j, i] * x[i] for i in 1:nx) +
                       sum(D[j, k] * input.u[k] for k in 1:nu) for j in 1:ny]...,
    ]
    compose(ODESystem(eqs, t, vcat(x...), [], name = name), [input, output])
end

StateSpace(A, B, C, D = nothing; kwargs...) = StateSpace(; A, B, C, D, kwargs...)
