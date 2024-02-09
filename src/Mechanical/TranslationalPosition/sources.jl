"""
    Force(; name)

Input signal acting as external force on a flange
"""
@mtkmodel Force begin
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(; use_support)
    @components begin
        f = RealInput(unit = u"N") # Accelerating force acting at flange (= -flange.tau)
    end
    @equations begin
        flange.f ~ -f.u
    end
end
