using ModelingToolkitStandardLibrary.Thermal, ModelingToolkit, OrdinaryDiffEq, Test
@parameters t

# Modelica example
begin
    @named mass1 = HeatCapacitor(C=15, T_start=373.15)
    @named mass2 = HeatCapacitor(C=15, T_start=273.15)
    @named conduction = ThermalConductor(G=10)
    @named Tsensor1 = TemperatureSensor() 
    @named Tsensor2 = TemperatureSensor()

    connections = [
        connect(mass1.port, conduction.port_a),
        connect(conduction.port_b, mass2.port),
        connect(mass1.port, Tsensor1.port),
        connect(mass2.port, Tsensor2.port),
    ]

    @named model = ODESystem(connections, t, systems=[mass1, mass2, conduction, Tsensor1, Tsensor2])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0, 3.0))
    sol = solve(prob, Rodas4())
end