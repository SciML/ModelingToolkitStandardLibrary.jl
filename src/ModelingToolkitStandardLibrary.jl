module ModelingToolkitStandardLibrary

"""
  @symcheck J > 0 || throw(ArgumentError("Expected `J` to be positive"))
  
Omits the check expression if the argument `J` is symbolic.
"""
macro symcheck(ex)
    ex.args[1].head === :call ||
        error("Expected an expresion on the form sym > val || error()")
    sym = ex.args[1].args[2]
    quote
        _issymbolic($(esc(sym))) || ($(esc(ex)))
    end
end

_issymbolic(x) = !(Symbolics.unwrap(x) isa Real)

include("Blocks/Blocks.jl")
include("Mechanical/Mechanical.jl")
include("Thermal/Thermal.jl")
include("Electrical/Electrical.jl")
include("Magnetic/Magnetic.jl")
include("Hydraulic/Hydraulic.jl")

end
