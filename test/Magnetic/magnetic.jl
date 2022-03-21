using ModelingToolkitStandardLibrary.Magnetic, ModelingToolkit, OrdinaryDiffEq, Test

@parameters t
@named ground = Ground()

@info "Testing basic magnetic components..."

