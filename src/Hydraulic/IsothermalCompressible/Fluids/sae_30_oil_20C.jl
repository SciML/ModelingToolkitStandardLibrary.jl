@named sae_30_oil_20C = FluidType()
const sae_30_oil_20C_type = typeof(get_fluid(sae_30_oil_20C))

# https://wiki.anton-paar.com/us-en/engine-oil/ at 20C
density(::sae_30_oil_20C_type) = 881.5
bulk_modulus(::sae_30_oil_20C_type) = 1.5e9
viscosity(::sae_30_oil_20C_type) = 0.23939