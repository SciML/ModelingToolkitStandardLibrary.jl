"""
    Force(; name)

Input signal acting as external force on a flange
"""
@mtkmodel Force begin
    @parameters begin
        s = 0
    end
    @components begin
        flange = Flange(; s = s)
        f = RealInput()
    end
    @equations begin
        flange.f ~ -f.u
    end
end
