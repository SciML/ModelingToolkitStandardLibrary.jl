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

@enum FluidType water sae_30_oil

# Base.float(x::FluidType) = Base.float(Int(x))

get_fluid(x::Float64) = get_fluid(Int(x))
get_fluid(x::Int) = get_fluid(FluidType(x))
get_fluid(x::FluidType) = Val(x)


density(fluid) = density(get_fluid(fluid))
@register_symbolic density(fluid)

bulk_modulus(fluid) = bulk_modulus(Val(fluid))
@register_symbolic bulk_modulus(fluid)

viscosity(fluid) = viscosity(get_fluid(fluid))
@register_symbolic viscosity(fluid)

density(::Val{water}) = 997 # [kg/m^3]
bulk_modulus(::Val{water}) = 2.09e9 # [Pa]
viscosity(::Val{water}) = 0.0010016 # [Pa*s] or [kg/m-s]

# https://wiki.anton-paar.com/us-en/engine-oil/ at 20C
density(::Val{sae_30_oil}) = 881.5
bulk_modulus(::Val{sae_30_oil}) = 1.5e9
viscosity(::Val{sae_30_oil}) = 0.23939


density(fluid, p) = density(fluid)*(1 + p/bulk_modulus(fluid))


@enum ShapeType circle triangle square rectangle

# Base.float(x::ShapeType) = Base.float(Int(x))

get_shape(x::Float64) = get_shape(Int(x))
get_shape(x::Int) = get_shape(ShapeType(x))
get_shape(x::FluidType) = Val(x)

friction_factor(shape) = friction_factor(get_shape(shape))
@register_symbolic friction_factor(shape)

friction_factor(::Val{circle}) = 64
