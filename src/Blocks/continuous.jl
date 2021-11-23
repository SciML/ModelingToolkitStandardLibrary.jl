# TODO: remove initial values for all inputs once IO handling in MTK is in place
"""
    Constant(val; name)

Outputs a constant value `val`.
"""
function Constant(val; name)
    @variables y(t)=val [output=true]
    @parameters val=val
    eqs = [
        y ~ val
    ]
    ODESystem(eqs, t, name=name)
end

"""
    Integrator(; k=1, name)

Outputs `y = ∫k*u dt`, corresponding to the transfer function `1/s`.
"""
function Integrator(; k=1, name)
    @variables x(t)=0 u(t)=0 [input=true] y(t)=0 [output=true]
    @parameters k=k
    eqs = [
        Dₜ(x) ~ k*u
        y ~ x
    ]
    ODESystem(eqs, t, name=name)
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
function Derivative(; k=1, T, name)
    @variables x(t)=0 u(t)=0 [input=true] y(t)=0 [output=true]
    @parameters T=T k=k
    eqs = [
        Dₜ(x) ~ (u - x) / T
        y ~ (k/T)*(u - x)
    ]
    ODESystem(eqs, t, name=name)
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
    @variables x(t)=0 u(t)=0 [input=true] y(t) [output=true]
    @parameters T=T k=k
    eqs = [
        Dₜ(x) ~ (-x + k*u) / T
        y ~ x
    ]
    ODESystem(eqs, t, name=name)
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
    @variables x(t)=0 xd(t)=0 u(t)=0 [input=true] y(t) [output=true]
    @parameters k=k w=w d=d
    eqs = [
        Dₜ(x) ~ xd
        Dₜ(xd) ~ w*(w*(k*u - x) - 2*d*xd)
        y ~ x
    ]
    ODESystem(eqs, t, name=name)
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
`gains`: If `gains = true`, `Ti` and `Td` will be interpreted as gains with a fundamental PID transfer function on parallel form `ki=Ti, kd=Td, k + ki/s + kd*s`
"""
function PID(; k, Ti=false, Td=false, wp=1, wd=1,
    Ni    = Ti == 0 ? Inf : √(max(Td / Ti, 1e-6)),
    Nd    = 12,
    y_max = Inf,
    y_min = y_max > 0 ? -y_max : -Inf,
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
    y_max ≥ y_min || throw(ArgumentError("y_min must be smaller than y_max"))

    @variables x(t)=0 u_r(t)=0 [input=true] u_y(t)=0 [input=true] y(t) [output=true] e(t)=0 ep(t)=0 ed(t)=0 ea(t)=0
    
    
    @named D = Derivative(k = Td, T = Td/Nd) # NOTE: consider T = max(Td/Nd, 100eps()), but currently errors since a symbolic variable appears in a boolean expression in `max`.
    if isequal(Ti, false)
        @named I = Gain(false)
    else
        @named I = Integrator(k = 1/Ti)
    end
    @named sat = Saturation(; y_min, y_max)
    derivative_action = Td > 0
    @parameters k=k Td=Td wp=wp wd=wd Ni=Ni Nd=Nd # TODO: move this line above the subsystem definitions when symbolic default values for parameters works. https://github.com/SciML/ModelingToolkit.jl/issues/1013
    # NOTE: Ti is not included as a parameter since we cannot support setting it to false after this constructor is called. Maybe Integrator can be tested with Ti = false setting k to 0 with IfElse?
    
    eqs = [
        e ~ u_r - u_y # Control error
        ep ~ wp*u_r - u_y  # Control error for proportional part with setpoint weight
        ea ~ sat.y - sat.u # Actuator error due to saturation
        I.u ~ e + 1/(k*Ni)*ea  # Connect integrator block. The integrator integrates the control error and the anti-wind up tracking. Note the apparent tracking time constant 1/(k*Ni), since k appears after the integration and 1/Ti appears in the integrator block, the final tracking gain will be 1/(Ti*Ni) 
        sat.u ~ derivative_action ? k*(ep + I.y + D.y) : k*(ep + I.y) # unsaturated output = P + I + D
        y ~ sat.y
    ]
    systems = [I, sat]
    if derivative_action
        push!(eqs, ed ~ wd*u_r - u_y)
        push!(eqs, D.u ~ ed) # Connect derivative block
        push!(systems, D)
    end
    ODESystem(eqs, t, name=name, systems=systems)
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
function StateSpace(A, B, C, D=0; x0=zeros(size(A,1)), name)
    nx = size(A,1)
    nu = size(B,2)
    ny = size(C,1)
    if nx == 0
        length(C) == length(B) == 0 || throw(ArgumentError("Dimension mismatch between A,B,C matrices"))
        return Gain(D; name=name)
    end
    if B isa AbstractVector
        B = reshape(B, length(B), 1)
    end
    if D == 0
        D = zeros(ny, nu)
    end
    @variables x[1:nx](t)=x0 u[1:nu](t)=0 [input=true] y[1:ny](t)=C*x0 [output=true]
    x = collect(x) # https://github.com/JuliaSymbolics/Symbolics.jl/issues/379
    u = collect(u)
    y = collect(y)
    # @parameters A=A B=B C=C D=D # This is buggy
    eqs = [
        Dₜ.(x) .~ A*x .+ B*u
        y      .~ C*x .+ D*u
    ]
    ODESystem(eqs, t, name=name)
end
