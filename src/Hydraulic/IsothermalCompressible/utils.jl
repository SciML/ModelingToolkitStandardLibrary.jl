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
- `p`: [Pa] pressure
- `dm`: [kg/s] mass flow
""" HydraulicPort


struct FluidType
    name::Symbol
    value::Int
end
FluidType(value; name) = FluidType(name, value)


Base.float(x::FluidType) = Base.float(x.value)
Base.eltype(x::FluidType) = Int
Base.promote_type(x::Type{FluidType}, y::Type{Int64}) = Int64
Base.promote_type(x::Type{FluidType}, y::Type{Float64}) = Float64
Base.convert(T::Type{Float64}, x::FluidType) = float(x)

Base.show(io::IO, x::FluidType) = show(io, MIME"text/plain"(), x)
function Base.show(io::IO, ::MIME"text/plain", x::FluidType)
    print(io, "$(x.name)::FluidType")
end

get_fluid(x::Float64) = Val(round(Int, x))
get_fluid(x::Int) = Val(x)
get_fluid(x::FluidType) = Val(x.value)


module Fluids
import ..IsothermalCompressible: FluidType
using ModelingToolkit
@named water = FluidType(0)
@named sae_30_oil = FluidType(1)
end



density(fluid) = density(get_fluid(fluid))
@register_symbolic density(fluid)

bulk_modulus(fluid) = bulk_modulus(get_fluid(fluid))
@register_symbolic bulk_modulus(fluid)

viscosity(fluid) = viscosity(get_fluid(fluid))
@register_symbolic viscosity(fluid)


density(::Val{Fluids.water.value}) = 997 # [kg/m^3]
bulk_modulus(::Val{Fluids.water.value}) = 2.09e9 # [Pa]
viscosity(::Val{Fluids.water.value}) = 0.0010016 # [Pa*s] or [kg/m-s]

# https://wiki.anton-paar.com/us-en/engine-oil/ at 20C
density(::Val{Fluids.sae_30_oil.value}) = 881.5
bulk_modulus(::Val{Fluids.sae_30_oil.value}) = 1.5e9
viscosity(::Val{Fluids.sae_30_oil.value}) = 0.23939


density(fluid, p) = density(fluid)*(1 + p/bulk_modulus(fluid))




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
@named triangle = ShapeType(1)
end


friction_factor(shape) = friction_factor(get_shape(shape))
@register_symbolic friction_factor(shape)

friction_factor(::Val{Shapes.circle.value}) = 64
