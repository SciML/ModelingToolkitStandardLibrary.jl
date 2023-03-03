@connector function HydraulicPort(; name, p_int)
    pars = []
    vars = @variables begin
        p(t) = p_int
        dm(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, pars, name = name, defaults = [dm => 0])
end
Base.@doc """
    HydraulicPort(;name)

port for hydraulic components

# States:
- `p`: [Pa] gage pressure
- `dm`: [kg/s] mass flow
""" HydraulicPort


struct FluidType
    name::Symbol
    value::Float64
end
FluidType(; name) = FluidType(name, float(hash(name)))


Base.float(x::FluidType) = Base.float(x.value)
Base.eltype(x::FluidType) = Float64
# Base.promote_type(x::Type{FluidType}, y::Type{Int64}) = Int64
Base.promote_type(x::Type{FluidType}, y::Type{Float64}) = Float64
Base.convert(T::Type{Float64}, x::FluidType) = float(x)

Base.show(io::IO, x::FluidType) = show(io, MIME"text/plain"(), x)
function Base.show(io::IO, ::MIME"text/plain", x::FluidType)
    print(io, "$(x.name)::FluidType")
end

get_fluid(x::Float64) = Val(x)
get_fluid(x::FluidType) = Val(x.value)


density(fluid) = density(get_fluid(fluid))
@register_symbolic density(fluid)

bulk_modulus(fluid) = bulk_modulus(get_fluid(fluid))
@register_symbolic bulk_modulus(fluid)

viscosity(fluid) = viscosity(get_fluid(fluid))
@register_symbolic viscosity(fluid)


density(fluid, p) = density(fluid)*(1 + p/bulk_modulus(fluid))

# Fluids
include("Fluids/water_20C.jl")
include("Fluids/sae_30_oil_20C.jl")




struct ShapeType
    name::Symbol
    value::Int
end
ShapeType(value; name) = ShapeType(name, value)


Base.float(x::ShapeType) = Base.float(x.value)
Base.eltype(x::ShapeType) = Int
Base.promote_type(x::Type{ShapeType}, y::Type{Int64}) = Int64
Base.promote_type(x::Type{ShapeType}, y::Type{Float64}) = Float64
Base.convert(T::Type{Float64}, x::ShapeType) = float(x)

Base.show(io::IO, x::ShapeType) = show(io, MIME"text/plain"(), x)
function Base.show(io::IO, ::MIME"text/plain", x::ShapeType)
    print(io, "$(x.name)::ShapeType")
end

get_shape(x::Float64) = Val(round(Int, x))
get_shape(x::Int) = Val(x)
get_shape(x::ShapeType) = Val(x.value)


module Shapes
import ..IsothermalCompressible: ShapeType
using ModelingToolkit
@named circle = ShapeType(0)
@named square = ShapeType(1)
@named triangle = ShapeType(2)
@named rectangle_r2 = ShapeType(3)
@named rectangle_r3 = ShapeType(4)
@named rectangle_r4 = ShapeType(5)
@named rectangle_r8 = ShapeType(6)
@named rectangle_rinf = ShapeType(7)
end

shape_factor(shape) = shape_factor(get_shape(shape))
@register_symbolic shape_factor(shape)

shape_factor(::Val{Shapes.circle.value}) = 64
shape_factor(::Val{Shapes.square.value}) = 57
shape_factor(::Val{Shapes.triangle.value}) = 53
shape_factor(::Val{Shapes.rectangle_r2.value}) = 62
shape_factor(::Val{Shapes.rectangle_r3.value}) = 69
shape_factor(::Val{Shapes.rectangle_r4.value}) = 73
shape_factor(::Val{Shapes.rectangle_r8.value}) = 82
shape_factor(::Val{Shapes.rectangle_rinf.value}) = 96



function friction_factor(dm, area, d_h, rho, mu, shape_factor)

    u = abs(dm)/(rho*area)
    Re = maximum([rho*u*d_h/mu, 1]) 

    f_laminar = shape_factor/Re


    f_turbulent = (0.79*log(Re)-1.64)^(-2)

    f = transition(2000, 3000, f_laminar, f_turbulent, Re)

    return f
end
@register_symbolic friction_factor(dm, area, d_h, rho, mu, shape_factor)



function transition(x1,x2,y1,y2,x)

    if x < x1
        return y1
    elseif x > x2
        return y2
    else
        u = (x-x1)/(x2-x1)
        blend = 3*u^2 - 2*u^3;
        return (1-blend)*y1 + blend*y2;
    end
    
end


