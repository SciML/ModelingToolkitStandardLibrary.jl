module Machines

using ModelingToolkitBase, Symbolics
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ...Electrical: Pin
using ...Mechanical.Rotational: Flange, Support

export PMSM
include("pmsm.jl")

end #module
