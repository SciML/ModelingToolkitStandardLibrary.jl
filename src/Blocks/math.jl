import Symbolics.scalarize

"""
    Gain(k; name)

Outputs `y = k*u`. `k` can be a scalar or an array.
"""
function Gain(k=1; name)
    @variables u(t)=0 [input=true] y(t) [output=true]
    @parameters k=k
    eqs = [
        y ~ k*u
    ]
    ODESystem(eqs, t, name=name)
end

function Gain(K::AbstractArray; name)
    ny,nu = size(K, 1), size(K, 2)
    @variables u[1:nu](t)=0 [input=true] y[1:ny](t)=0 [output=true]
    eqs = y .~ K*u
    ODESystem(eqs, t, name=name)
end

"""
    Sum(n::Int; name)
    Sum(k::AbstractVector; name)

Creates a summing block that sums `n` inputs, `y = sum(u[i] for i ∈ 1:n)`.
A vector of summing coefficients `k` can also be provided, i.e., `y = sum(k[i]u[i] for i ∈ 1:n)`.
A block that subtracts one signal from another can thus be created by `@named sub = Sum([1, -1])`.
"""
function Sum(n::Int; name)
    @variables u[1:n](t)=0 [input=true] y(t)=0 [output=true]
    eqs = [y ~ scalarize(sum(u))]
    ODESystem(eqs, t, name=name)
end

function Sum(k::AbstractVector; name)
    n = length(k)
    @variables u[1:n](t)=0 [input=true] y(t)=0 [output=true]
    eqs = [y ~ scalarize(sum(k[i]*u[i] for i ∈ 1:n))]
    ODESystem(eqs, t, name=name)
end

function Product(n::Int=2; name)
    @variables u[1:n](t)=0 [input=true] y(t)=0 [output=true]
    eqs = [y ~ scalarize(prod(u))]
    ODESystem(eqs, t, name=name)
end