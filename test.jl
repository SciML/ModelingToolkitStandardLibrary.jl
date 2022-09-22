using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkit
using OrdinaryDiffEq

@named link = Link(; m = 1, l = 1, I = 1, g = -10)

eqs = [link.TX1.v ~ 0
       link.TY1.v ~ 0
       link.TX1.f ~ 0
       link.TY1.f ~ 0]

@named model = ODESystem(eqs, t, [], []; systems = [link])

sys = structural_simplify(model)
