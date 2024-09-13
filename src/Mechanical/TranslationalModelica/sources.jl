"""
    Force(; name, use_support = false)

Input signal acting as external force on a flange
"""
@mtkmodel Force begin
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = false)
    @components begin
        f = RealInput() # Accelerating force acting at flange (= -flange.tau)
    end
    @equations begin
        flange.f ~ -f.u
    end
end

@mtkmodel Position begin
    @extend (s,) = ptf = PartialElementaryOneFlangeAndSupport2()
    @structural_parameters begin
        exact = false
    end
    @parameters begin
        f_crit = 50
    end
    @variables begin
        v(t)
        a(t)
    end
    @components begin
        s_ref = RealInput()
    end
    begin
        w_crit = 2π * f_crit
        af = 1.3617
        bf = 0.6180
    end
    @equations begin
        if exact
            s ~ s_ref.u
        else
            a ~ ((s_ref.u - s) * w_crit - af * v) * (w_crit / bf)
        end
        v ~ D(s)
        a ~ D(v)
    end
end
