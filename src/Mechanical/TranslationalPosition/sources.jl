"""
    Force(; name, use_support = false)

Input signal acting as external force on a flange
"""
@mtkmodel Force begin
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = false)
    @parameters begin
        use_support = false
        s = 0
    end
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = use_support isa Bool ? use_support :
                      ModelingToolkit.getdefault(use_support))
    @components begin
        flange = Flange(; s = s)
        f = RealInput()
    end
    @equations begin
        flange.f ~ -f.u
    end
end
