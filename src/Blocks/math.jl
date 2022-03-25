"""
    Gain(k; name)

Outputs `y = k*u`. `k` can be a scalar or an array.
"""
function Gain(k=1; name)
    @named u = RealInput()
    @named y = RealOutput()
    pars = @parameters k=k
    eqs = [
        y.u ~ k*u.u
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [y, u])
end

function Gain(K::AbstractArray; name) # FIXME:
    ny,nu = size(K, 1), size(K, 2)
    @variables u[1:nu](t)=0 [input=true] y[1:ny](t)=0 [output=true]
    u = collect(u)
    y = collect(y)
    eqs = y .~ K*u
    ODESystem(eqs, t, [u, y], [], name=name)
end

"""
    Sum(n::Int; name)
    Sum(k::AbstractVector; name)

Creates a summing block that sums `n` inputs, `y = sum(u[i] for i ∈ 1:n)`.
A vector of summing coefficients `k` can also be provided, i.e., `y = sum(k[i]u[i] for i ∈ 1:n)`.
A block that subtracts one signal from another can thus be created by `@named sub = Sum([1, -1])`.
"""
function Sum(n::Int; name) # FIXME:
    @variables u[1:n](t)=0 [input=true] y(t)=0 [output=true]
    u = collect(u)
    eqs = [y ~ sum(u)]
    ODESystem(eqs, t, [u, y], name=name)
end

function Sum(k::AbstractVector; name)# FIXME:
    n = length(k)
    @variables u[1:n](t)=0 [input=true] y(t)=0 [output=true]
    u = collect(u)
    eqs = [y ~ sum(k[i]*u[i] for i ∈ 1:n)]
    ODESystem(eqs, t, [u, y], name=name)
end

function Product(n::Int=2; name) # FIXME:
    @variables u[1:n](t)=0 [input=true] y(t)=0 [output=true]
    u = collect(u)
    eqs = [y ~ prod(u)]
    ODESystem(eqs, t, [u, y], name=name)
end