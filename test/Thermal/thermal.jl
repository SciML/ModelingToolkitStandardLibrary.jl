using ModelingToolkitStandardLibrary.Thermal, ModelingToolkit, OrdinaryDiffEq, Test
@parameters t

# Test HeatCapacitor, TemperatureSensor, RelativeTemperatureSensor, FixedTemperature
@testset "Heat systems" begin
    T, C, G = 10.0, 10.0, 10.0
    @variables final_T(t)
    @named mass1         = HeatCapacitor(C=C)
    @named mass2         = HeatCapacitor(C=C)
    @named th_conductor  = ThermalConductor(G=G)
    @named reltem_sensor = RelativeTemperatureSensor()
    @named T_sensor1     = TemperatureSensor()
    @named T_sensor2     = TemperatureSensor()
    @named tem_src       = FixedTemperature(T=T)

    @info "Building a single-body system..."
    eqs = [
        connect(mass1.a, th_conductor.a)
        connect(th_conductor.b, reltem_sensor.a)
        connect(reltem_sensor.b, tem_src.b)
    ]
    @named h1 = ODESystem(eqs, t, systems=[mass1, reltem_sensor, tem_src, th_conductor])
    sys = structural_simplify(h1)

    u0 = [
        mass1.T => 2.0
    ]
    prob = ODEProblem(sys, u0, (0, 2.0))
    sol = solve(prob, Rodas4())
    
    # Check if Relative temperature sensor reads the temperature of heat capacitor
    # when connected to a thermal conductor and a fixed temperature source 
    @test sol[reltem_sensor.T] + sol[tem_src.b.T] == sol[mass1.T] + sol[th_conductor.dT] 
    
    @info "Building a two-body system..."
    eqs = [
        connect(T_sensor1.a, mass1.a, th_conductor.a)
        connect(th_conductor.b, mass2.a, T_sensor2.a)
        final_T ~ (mass1.C * mass1.T + mass2.C * mass2.T) / 
        (mass1.C + mass2.C)
    ]
    @named h2 = ODESystem(eqs, t, [final_T], [],
                        systems=[mass1, mass2, T_sensor1, T_sensor2, th_conductor])
    sys = structural_simplify(h2)

    u0 = [
        mass1.T => 1.0
        mass2.T => 10.0
        final_T => 12
    ]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())
    
    m1, m2 = sol.u[end]
    @test m1 ≈ m2 atol=1e-1
    mass_T = reduce(hcat, sol.u)
    @test sol[T_sensor1.T] == mass_T[1, :]
    @test sol[T_sensor2.T] == mass_T[2, :]
end

# Test HeatFlowSensor, FixedHeatFlow, ThermalResistor, ThermalConductor, ThermalGround
@testset "Heat flow system" begin
    C, G, R = 10, 10, 10
    @named flow_src     = FixedHeatFlow(Q_flow=50, alpha=100)
    @named mass1        = HeatCapacitor(C=C)
    @named hf_sensor1   = HeatFlowSensor()
    @named hf_sensor2   = HeatFlowSensor()
    @named th_conductor = ThermalConductor(G=G)
    @named th_resistor  = ThermalResistor(R=R) 
    @named th_ground    = ThermalGround()

    @info "Building a heat-flow system..."
    eqs = [
        connect(mass1.a, th_resistor.a, th_conductor.a)
        connect(th_conductor.b, flow_src.b, hf_sensor1.a, hf_sensor2.a)
        connect(th_resistor.b, hf_sensor1.b, hf_sensor2.b, th_ground.a)
    ]
    @named h2 = ODESystem(eqs, t, 
                          systems=[mass1, hf_sensor1, hf_sensor2, 
                                   th_resistor, flow_src, th_ground, th_conductor])
    sys = structural_simplify(h2)

    u0 = [
        mass1.T            => 10.0
        th_resistor.Q_flow => 1.0
    ]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol  = solve(prob, Rodas4())
    
    @test sol[th_conductor.dT] .* G == sol[th_conductor.Q_flow]
    @test sol[th_conductor.Q_flow] ≈ sol[hf_sensor1.Q_flow] + sol[flow_src.b.Q_flow]

    @test sol[mass1.T] == sol[th_resistor.a.T]
    @test sol[th_resistor.dT] ./ R ≈ sol[th_resistor.Q_flow]

end

# Tests ConvectiveResistor and includes FixedTemperature and ThermalResistor
@testset "Piston cylinder wall" begin
    Tᵧ, Tᵪ = 1000, 10 # ᵧ -> gas and ᵪ -> coolant
    Rᵧ, Rᵪ = 50e-4, 10e-4 # R = 1/h; h is convection co-efficient
    R_wall = 1.5e-4
    @named coolant     = ConvectiveResistor(R=Rᵪ)
    @named gas         = ConvectiveResistor(R=Rᵧ)
    @named wall        = ThermalResistor(R=R_wall)
    @named gas_tem     = FixedTemperature(T=Tᵧ)
    @named coolant_tem = FixedTemperature(T=Tᵪ)

    @info "Building a piston-cylinder..."
    eqs = [
        connect(gas_tem.b, gas.solidport)
        connect(gas.fluidport, wall.a)
        connect(wall.b, coolant.fluidport)
        connect(coolant.solidport, coolant_tem.b)
    ]
    @named piston = ODESystem(eqs, t, systems=[gas_tem, wall, gas, coolant, coolant_tem])
    sys = structural_simplify(piston)

    u0 = [
        coolant.dT  => 5.0
        wall.Q_flow => 10.0
    ]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())
    
    # Heat-flow-rate is equal in magnitude 
    # and opposite in direction
    @test sol[gas.Q_flow] + sol[coolant.Q_flow] == zeros(length(sol))
end

# Test ConvectiveConductor, BodyRadiation
@testset "Radiator system" begin
    Tᵧ, Tᵪ = 1000, 10 # ᵧ -> gas and ᵪ -> coolant
    R_wall = 10
    G = 0.04
    σ   = 5.6703744191844294e-8 # Stefan-Boltzmann constant

    @named base        = ThermalResistor(R=R_wall)
    @named gas_tem     = FixedTemperature(T=Tᵧ)
    @named coolant_tem = FixedTemperature(T=Tᵪ)
    @named radiator    = BodyRadiation(G=G)
    @named ground      = ThermalGround()
    @named dissipator  = ConvectiveConductor(G=10)
    
    @info "Building a radiator..."
    eqs = [
        connect(gas_tem.b, radiator.a, base.a, dissipator.solidport)
        connect(base.b, radiator.b, coolant_tem.b, dissipator.fluidport)
    ]
    @named rad = ODESystem(eqs, t, systems=[base, gas_tem, radiator, dissipator, coolant_tem])
    sys = structural_simplify(rad)

    u0 = [
        base.Q_flow => 10
        dissipator.Q_flow => 10
    ]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())

    @test sol[dissipator.dT] == sol[radiator.a.T] - sol[radiator.b.T]
    rad_Q_flow = G*σ*(Tᵧ^4 - Tᵪ^4)
    @test sol[radiator.Q_flow] == fill(rad_Q_flow, length(sol[radiator.Q_flow]))
end

@testset "Thermal Collector" begin
    @named flow_src    = FixedHeatFlow(Q_flow=50, alpha=100)
    @named hf_sensor   = HeatFlowSensor()
    @named th_ground   = ThermalGround()
    @named collector   = ThermalCollector(N=2)
    @named th_resistor = ThermalResistor(R=10) 
    @named tem_src     = FixedTemperature(T=10)

    @info "Building a heat collector..."
    eqs = [
        connect(flow_src.b, collector.hp1, th_resistor.a)
        connect(tem_src.b, collector.hp2)
        connect(hf_sensor.a, collector.collector_port)
        connect(hf_sensor.b, th_ground.a, th_resistor.b)
        ]
    @named coll = ODESystem(eqs, t, 
                            systems=[hf_sensor,flow_src, tem_src, 
                            collector, th_resistor])
    sys = structural_simplify(coll)

    u0 = [
        th_resistor.Q_flow => 1.0
    ]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol  = solve(prob, Rodas4())

    @test sol[collector.collector_port.Q_flow] + sol[collector.hp1.Q_flow] + sol[collector.hp2.Q_flow] ==
        zeros(length(sol[collector.collector_port.Q_flow]))
    @test sol[collector.collector_port.T] == sol[collector.hp1.T] == sol[collector.hp2.T]
end