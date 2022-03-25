@connector function MagneticPort(;name)
    sts = @variables begin
        V_m(t) # [Wb] Magnetic potential at the port
        Phi(t), [connect=Flow] # [A] Magnetic flux flowing into the port"
    end
    ODESystem(Equation[], t, sts, []; name=name)
end
Base.@doc "Port for a Magnetic system." MagneticPort