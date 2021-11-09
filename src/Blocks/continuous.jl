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
