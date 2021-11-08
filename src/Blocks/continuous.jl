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
