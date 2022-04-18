"""
    Integrator(; k=1, name)

Outputs `y = ∫k*u dt`, corresponding to the transfer function `1/s`.
"""
function Integrator(;name, k=1, x0=0.0)
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t)=x0
    pars = @parameters k=k
    eqs = [
        D(x) ~ k * u
        y ~ x
    ]
    extend(ODESystem(eqs, t, sts, pars; name=name), siso)
end

"""
    Derivative(; k=1, T, name)

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
"""
function Derivative(; name, k=1, T=10, x0=0)
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t)=x0
    pars = @parameters T=T k=k
    eqs = [
        D(x) ~ (u - x) / T
        y ~ (k / T) * (u - x)
    ]
    extend(ODESystem(eqs, t, sts, pars; name=name), siso)
end

"""
    FirstOrder(; k=1, T, name)

A first-order filter with a single real pole in `s = -T` and gain `k`. The transfer function
is given by `Y(s)/U(s) = `
```
   k   
───────
sT + 1
```
"""
function FirstOrder(; k=1, T, name)
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t)=0
    pars = @parameters T=T k=k
    eqs = [
        D(x) ~ (u - x) / T
        y ~ x
    ]
    extend(ODESystem(eqs, t, sts, pars; name=name), siso)
end

"""
    SecondOrder(; k=1, w, d, name)

A second-order filter with gain `k`, a bandwidth of `w` rad/s and relative damping `d`. The transfer function
is given by `Y(s)/U(s) = `
```
      k*w^2   
─────────────────
s² + 2d*w*s + w^2
```
Critical damping corresponds to `d=1`, which yields the fastest step response without overshoot, d < 1` results in an under-damped filter while `d > 1` results in an over-damped filter.
`d = 1/√2` corresponds to a Butterworth filter of order 2 (maximally flat frequency response).
"""
function SecondOrder(; k=1, w, d, name)
    @named siso = SISO()
    @unpack u, y = siso
    sts = @variables x(t)=0 xd(t)=0
    pars = @parameters k=k w=w d=d
    eqs = [
        D(x) ~ xd
        D(xd) ~ w*(w*(k*u - x) - 2*d*xd)
        y ~ x
    ]
    extend(ODESystem(eqs, t, sts, pars; name=name), siso)
end

"""
PI-controller without actuator saturation and anti-windup measure.
"""
function PI(;name, k=1, T=1, x_start=0)
    @named e = RealInput() # control error
    @named u = RealOutput() # control signal
    @variables x(t)=x_start
    T > 0 || error("Time constant `T` has to be strictly positive")
    pars = @parameters k=k T=T
    eqs = [
        D(x) ~ 1 / T * e.u
        u.u ~ k * (x + e.u)
    ]
    compose(ODESystem(eqs, t, [x], pars; name=name), [e, u])
end

"""
Text-book version of a PID-controller.
"""
function PID(;name, k=1, Ti=1, Td=1, Nd=10, xi_start=0, xd_start=0)
    @named err_input = RealInput() # control error
    @named ctr_output = RealOutput() # control signal
    Ti > 0 || error("Time constant `Ti` has to be strictly positive")
    Td > 0 || error("Time constant `Td` has to be strictly positive")
    Nd > 0 || error("`Nd` has to be strictly positive")
    @named gain = Gain(k)
    @named int = Integrator(k=1/Ti, x0=xi_start)
    @named der = Derivative(k=1/Td, T=1/Nd, x0=xd_start)
    @named add = Add3()
    eqs = [
        connect(err_input, add.input1),
        connect(err_input, int.input),
        connect(err_input, der.input),
        connect(int.output, add.input2),
        connect(der.output, add.input3),
        connect(add.output, gain.input),
        connect(gain.output, ctr_output)
    ]
    ODESystem(eqs, t, [], []; name=name, systems=[gain, int, der, add, err_input, ctr_output])
end

"""
PI-controller with actuator saturation and anti-windup measure.
"""
function LimPI(;name, k=1, T=1, u_max=1, u_min=-u_max, Ta=1)
    @named e = RealInput() # control error
    @named u = RealOutput() # control signal
    @variables x(t)=0.0 u_star(t)=0.0
    Ta > 0 || error("Time constant `Ta` has to be strictly positive")
    T > 0 || error("Time constant `T` has to be strictly positive")
    pars = @parameters k=k T=T u_max=u_max u_min=u_min
    eqs = [
        D(x) ~ e.u * k / T + 1 / Ta * (-u_star + u.u)
        u.u ~ max(min(u_star, u_max), u_min)      
        u_star ~ x + k * e.u
    ]
    compose(ODESystem(eqs, t, [x, u_star], pars; name=name), [e, u])
end

"""
    PID(; k, Ti=false, Td=false, wp=1, wd=1, Ni, Nd=12, y_max=Inf, y_min=-y_max, gains = false, name)

Proportional-Integral-Derivative (PID) controller with output saturation, set-point weighting and integrator anti-windup.

This block has inputs
- `u_r` reference/set point
- `u_y` measurement signal
and output `y` corresponding to the control signal.

The equation for the control signal is roughly
```
k(ep + 1/Ti * ∫e + 1/Td * d/dt(ed))
e = u_r - u_y
ep = wp*u_r - u_y
ed = wd*u_r - u_y
```
where the transfer function for the derivative includes additional filtering, see `? Derivative` for more details.

- `k`: Proportional gain
- `Ti`: Integrator time constant. Set to `false` to turn off integral action.
- `Td`: Derivative time constant. Set to `false` to turn off derivative action.
- `wp`: Set-point weighting in the proportional part.
- `wd`: Set-point weighting in the derivative part.
- `Nd`: Derivative limit, limits the derivative gain to Nd/Td. Reasonable values are ∈ [8, 20]. A higher value gives a better approximation of an ideal derivative at the expense of higher noise amplification.
- `Ni`: `Ni*Ti` controls the time constant `Tₜ` of anti-windup tracking. A common (default) choice is `Tₜ = √(Ti*Td)` which is realized by `Ni = √(Td / Ti)`. Anti-windup can be effectively turned off by setting `Ni = Inf`.
` `gains`: If `gains = true`, `Ti` and `Td` will be interpreted as gains with a fundamental PID transfer function on parallel form `ki=Ti, kd=Td, k + ki/s + kd*s`
"""
function LimPID(; k, Ti=false, Td=false, wp=1, wd=1,
    Ni    = Ti == 0 ? Inf : √(max(Td / Ti, 1e-6)),
    Nd    = 12,
    u_max = Inf,
    u_min = u_max > 0 ? -u_max : -Inf,
    gains = false,
    name
    )
    if gains
        Ti = k / Ti
        Td = Td / k
    end
    0 ≤ wp ≤ 1 || throw(ArgumentError("wp out of bounds, got $(wp) but expected wp ∈ [0, 1]"))
    0 ≤ wd ≤ 1 || throw(ArgumentError("wd out of bounds, got $(wd) but expected wd ∈ [0, 1]"))
    Ti ≥ 0     || throw(ArgumentError("Ti out of bounds, got $(Ti) but expected Ti ≥ 0"))
    Td ≥ 0     || throw(ArgumentError("Td out of bounds, got $(Td) but expected Td ≥ 0"))
    u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))

    @named r = RealInput() # reference
    @named y = RealInput() # measurement
    @named u = RealOutput() # control signal

    sts = @variables x(t)=0 e(t)=0 ep(t)=0 ed(t)=0 ea(t)=0
   
    @named D = Derivative(k = Td, T = Td/Nd) # NOTE: consider T = max(Td/Nd, 100eps()), but currently errors since a symbolic variable appears in a boolean expression in `max`.
    if isequal(Ti, false)
        @named I = Gain(1)
    else
        @named I = Integrator(k = 1/Ti)
    end
    @named sat = Limiter(; y_min=y_min, y_max=y_max)
    derivative_action = Td > 0
    pars = @parameters k=k Td=Td wp=wp wd=wd Ni=Ni Nd=Nd # TODO: move this line above the subsystem definitions when symbolic default values for parameters works. https://github.com/SciML/ModelingToolkit.jl/issues/1013
    # NOTE: Ti is not included as a parameter since we cannot support setting it to false after this constructor is called. Maybe Integrator can be tested with Ti = false setting k to 0 with IfElse?
    
    eqs = [
        e ~ r.u - y.u # Control error
        ep ~ wp * r.u - y.u  # Control error for proportional part with setpoint weight
        ea ~ sat.y.u - sat.u.u # Actuator error due to saturation
        I.u ~ e + 1 / (k * Ni) * ea  # Connect integrator block. The integrator integrates the control error and the anti-wind up tracking. Note the apparent tracking time constant 1/(k*Ni), since k appears after the integration and 1/Ti appears in the integrator block, the final tracking gain will be 1/(Ti*Ni) 
        sat.u ~ derivative_action ? k * (ep + I.y + D.y) : k * (ep + I.y) # unsaturated output = P + I + D
        y ~ sat.y
    ]
    systems = [I, sat]
    if derivative_action
        push!(eqs, ed ~ wd*u_r - u_y)
        push!(eqs, D.u ~ ed) # Connect derivative block
        push!(systems, D)
    end
    ODESystem(eqs, t, sts, pars, name=name, systems=systems)
end

"""
    StateSpace(A, B, C, D=0; x0=zeros(size(A,1)), name)

A linear, time-invariant state-space system on the form.
```
ẋ = Ax + Bu
y = Cx + Du
```
Transfer functions can also be simulated by converting them to a StateSpace form.
"""
function StateSpace(;A, B, C, D=nothing, x0=zeros(size(A,1)), name)
    nx, nu, ny = size(A,1), size(B,2), size(C,1)
    size(A,2) == nx || error("`A` has to be a square matrix.")
    size(B,1) == nx || error("`B` has to be of dimension ($nx x $nu).")
    size(C,2) == nx || error("`C` has to be of dimension ($ny x $nx).")
    if B isa AbstractVector
        B = reshape(B, length(B), 1)
    end
    if isnothing(D)
        D = zeros(ny, nu)
    else
        size(D) == (ny,nu) || error("`D` has to be of dimension ($ny x $nu).")
    end
    @named input = RealInput(nin=nu)
    @named output = RealOutput(nout=ny)
    @variables x[1:nx](t)=x0
    # pars = @parameters A=A B=B C=C D=D # This is buggy
    eqs = [ # FIXME: if array equations work
        [Differential(t)(x[i]) ~ sum(A[i,k] * x[k] for k in 1:nx) + sum(B[i,j] * input.u[j] for j in 1:nu) for i in 1:nx]..., # cannot use D here
        [output.u[j] ~ sum(C[j,i] * x[i] for i in 1:nx) + sum(D[j,k] * input.u[k] for k in 1:nu) for j in 1:ny]...,
    ]
    compose(ODESystem(eqs, t, vcat(x...), [], name=name), [input, output])
end
