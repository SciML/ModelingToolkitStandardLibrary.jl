@named water_20C = FluidType()
const water_20C_type = typeof(get_fluid(water_20C)) 

density(::water_20C_type) = 997             # [kg/m^3]
bulk_modulus(::water_20C_type) = 2.09e9     # [Pa]
viscosity(::water_20C_type) = 0.0010016     # [Pa*s] or [kg/m-s]