using ModelingToolkitStandardLibrary.Thermal, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Blocks: Constant, Step
using OrdinaryDiffEq: ReturnCode.Success

@parameters t
D = Differential(t)
#=
# Test HeatCapacitor, TemperatureSensor, RelativeTemperatureSensor, FixedTemperature
@testset "Heat systems" begin
    T, C, G = 10.0, 10.0, 10.0
    @variables final_T(t)
    @named mass1 = HeatCapacitor(C = C)
    @named mass2 = HeatCapacitor(C = C)
    @named th_conductor = ThermalConductor(G = G)
    @named reltem_sensor = RelativeTemperatureSensor()
    @named T_sensor1 = TemperatureSensor()
    @named T_sensor2 = TemperatureSensor()
    @named tem_src = FixedTemperature(T = T)

    @info "Building a single-body system..."
    eqs = [connect(mass1.port, th_conductor.port_a)
        connect(th_conductor.port_b, reltem_sensor.port_a)
        connect(reltem_sensor.port_b, tem_src.port)]
    @named h1 = ODESystem(eqs, t, systems = [mass1, reltem_sensor, tem_src, th_conductor])
    sys = structural_simplify(h1)

    u0 = [mass1.T => 2.0
        mass1.der_T => 1.0]
    prob = ODEProblem(sys, u0, (0, 2.0))
    sol = solve(prob, Tsit5())

    # Check if Relative temperature sensor reads the temperature of heat capacitor
    # when connected to a thermal conductor and a fixed temperature source
    @test sol.retcode == Success
    @test sol[reltem_sensor.T] + sol[tem_src.port.T] == sol[mass1.T] + sol[th_conductor.dT]

    @info "Building a two-body system..."
    eqs = [connect(T_sensor1.port, mass1.port, th_conductor.port_a)
        connect(th_conductor.port_b, mass2.port, T_sensor2.port)
        final_T ~ (mass1.C * mass1.T + mass2.C * mass2.T) /
                  (mass1.C + mass2.C)]
    @named h2 = ODESystem(eqs, t, [final_T], [],
        systems = [mass1, mass2, T_sensor1, T_sensor2, th_conductor])
    sys = structural_simplify(h2)

    u0 = [mass1.T => 1.0
        mass2.T => 10.0
        final_T => 12
        mass1.der_T => 1.0
        mass2.der_T => 1.0]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Tsit5())

    @test sol.retcode == Success
    m1, m2 = sol.u[end]
    @test m1≈m2 atol=1e-1
    mass_T = reduce(hcat, sol.u)
    @test sol[T_sensor1.T] == mass_T[1, :]
    @test sol[T_sensor2.T] == mass_T[2, :]
end

# Test HeatFlowSensor, FixedHeatFlow, ThermalResistor, ThermalConductor
@testset "Heat flow system" begin
    C, G, R = 10, 10, 10
    @named flow_src = FixedHeatFlow(Q_flow = 50, alpha = 100)
    @named mass1 = HeatCapacitor(C = C)
    @named hf_sensor1 = HeatFlowSensor()
    @named hf_sensor2 = HeatFlowSensor()
    @named th_conductor = ThermalConductor(G = G)
    @named th_resistor = ThermalResistor(R = R)
    @named th_ground = FixedTemperature(T = 0)

    @info "Building a heat-flow system..."
    eqs = [connect(mass1.port, th_resistor.port_a, th_conductor.port_a)
        connect(th_conductor.port_b, flow_src.port, hf_sensor1.port_a,
        hf_sensor2.port_a)
        connect(th_resistor.port_b, hf_sensor1.port_b, hf_sensor2.port_b,
        th_ground.port)]
    @named h2 = ODESystem(eqs, t,
        systems = [mass1, hf_sensor1, hf_sensor2,
            th_resistor, flow_src, th_ground, th_conductor])
    sys = structural_simplify(h2)

    u0 = [mass1.T => 10.0
        th_resistor.Q_flow => 1.0
        mass1.der_T => 1.0]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Tsit5())

    @test sol.retcode == Success
    @test sol[th_conductor.dT] .* G == sol[th_conductor.Q_flow]
    @test sol[th_conductor.Q_flow] ≈ sol[hf_sensor1.Q_flow] + sol[flow_src.port.Q_flow]

    @test sol[mass1.T] == sol[th_resistor.port_a.T]
    @test sol[th_resistor.dT] ./ R ≈ sol[th_resistor.Q_flow]
end

# Tests ConvectiveResistor and includes FixedTemperature and ThermalResistor
@testset "Piston cylinder wall" begin
    Tᵧ, Tᵪ = 1000, 10 # ᵧ -> gas and ᵪ -> coolant
    Rᵧ, Rᵪ = 50e-4, 10e-4 # R = 1/h; h is convection co-efficient
    R_wall = 1.5e-4
    @named coolant = ConvectiveResistor(R = Rᵪ)
    @named gas = ConvectiveResistor(R = Rᵧ)
    @named wall = ThermalResistor(R = R_wall)
    @named gas_tem = FixedTemperature(T = Tᵧ)
    @named coolant_tem = FixedTemperature(T = Tᵪ)

    @info "Building a piston-cylinder..."
    eqs = [connect(gas_tem.port, gas.solid)
        connect(gas.fluid, wall.port_a)
        connect(wall.port_b, coolant.fluid)
        connect(coolant.solid, coolant_tem.port)]
    @named piston = ODESystem(eqs, t, systems = [gas_tem, wall, gas, coolant, coolant_tem])
    sys = structural_simplify(piston)

    u0 = [coolant.dT => 5.0
        wall.Q_flow => 10.0]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())

    # Heat-flow-rate is equal in magnitude
    # and opposite in direction
    @test sol.retcode == Success
    @test sol[gas.Q_flow] + sol[coolant.Q_flow] == zeros(length(sol))
end

# Test ConvectiveConductor, BodyRadiation
@testset "Radiator system" begin
    T_gas, T_coolant = 1000, 10
    R_wall = 10
    G = 0.04
    σ = 5.6703744191844294e-8 # Stefan-Boltzmann constant

    @named base = ThermalResistor(R = R_wall)
    @named gas_tem = FixedTemperature(T = T_gas)
    @named coolant_tem = FixedTemperature(T = T_coolant)
    @named radiator = BodyRadiation(G = G)
    @named dissipator = ConvectiveConductor(G = 10)
    @named mass = HeatCapacitor(C = 10)

    @info "Building a radiator..."
    eqs = [connect(gas_tem.port, radiator.port_a, base.port_a, dissipator.solid, mass.port)
        connect(coolant_tem.port, base.port_b, radiator.port_b, dissipator.fluid)]
    @named rad = ODESystem(eqs, t,
        systems = [
            base,
            gas_tem,
            radiator,
            dissipator,
            coolant_tem,
            mass,
        ])
    sys = structural_simplify(rad)

    u0 = [base.Q_flow => 10
        dissipator.Q_flow => 10
        mass.T => T_gas]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    @test sol[dissipator.dT] == sol[radiator.port_a.T] - sol[radiator.port_b.T]
    rad_Q_flow = G * σ * (T_gas^4 - T_coolant^4)
    @test sol[radiator.Q_flow] == fill(rad_Q_flow, length(sol[radiator.Q_flow]))
end

@testset "Thermal Collector" begin
    @named flow_src = FixedHeatFlow(Q_flow = 50, alpha = 100)
    @named hf_sensor = HeatFlowSensor()
    @named collector = ThermalCollector(m = 2)
    @named th_resistor = ThermalResistor(R = 10)
    @named tem_src = FixedTemperature(T = 10)
    @named mass = HeatCapacitor(C = 10)

    @info "Building a heat collector..."
    eqs = [connect(flow_src.port, collector.port_a1, th_resistor.port_a)
        connect(tem_src.port, collector.port_a2)
        connect(hf_sensor.port_a, collector.port_b)
        connect(hf_sensor.port_b, mass.port, th_resistor.port_b)]
    @named coll = ODESystem(eqs, t,
        systems = [hf_sensor, flow_src, tem_src,
            collector, th_resistor, mass])
    sys = structural_simplify(coll)

    u0 = [
        th_resistor.Q_flow => 1.0,
        mass.T => 0.0,
    ]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    @test sol[collector.port_b.Q_flow] + sol[collector.port_a1.Q_flow] +
          sol[collector.port_a2.Q_flow] ==
          zeros(length(sol[collector.port_b.Q_flow]))
    @test sol[collector.port_b.T] == sol[collector.port_a1.T] == sol[collector.port_a2.T]
end
=#
# https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Thermal.HeatTransfer.Examples.Motor.html
@testset "demo" begin
    k2c(T) = T - 273.15
    T_amb = 293.15
    @named windingLosses = PrescribedHeatFlow(T_ref = k2c(95), alpha = 3.03e-3)
    @named winding = HeatCapacitor(C = 2500, T = T_amb)
    @named T_winding = TemperatureSensor()
    @named winding2core = ThermalConductor(G = 10)
    @named coreLosses = PrescribedHeatFlow()
    @named core = HeatCapacitor(C = 25000, T = T_amb)
    @named T_core = TemperatureSensor()
    @named convection = ConvectiveConductor(G = 25)
    @named environment = PrescribedTemperature()
    @named amb = Constant(k = T_amb)
    @named core_losses_const = Constant(k = 500)
    @named winding_losses = Step(height = 900, offset = 100, start_time = 360,
        duration = Inf, smooth = false)
    connections = [connect(windingLosses.port, winding.port)
        connect(coreLosses.port, core.port)
        connect(winding.port, winding2core.port_a)
        connect(winding2core.port_b, core.port)
        connect(winding.port, T_winding.port)
        connect(core.port, T_core.port)
        connect(winding2core.port_b, convection.solid)
        connect(convection.fluid, environment.port)
        connect(amb.output, environment.T)
        connect(winding_losses.output, windingLosses.Q_flow)
        connect(core_losses_const.output, coreLosses.Q_flow)]

    @named model = ODESystem(connections, t,
        systems = [
            windingLosses, winding, T_winding, winding2core,
            coreLosses, core,
            T_core, convection, environment, amb, core_losses_const,
            winding_losses])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0, 720.0))
    sol = solve(prob, Rodas4())

    # plot(sol; vars=[T_winding.T, T_core.T])
    @test sol.retcode == Success
    @test sol[T_winding.T] == sol[winding.T]
    @test sol[T_core.T] == sol[core.T]
    @test sol[-core.port.Q_flow] ==
          sol[coreLosses.port.Q_flow + convection.solid.Q_flow + winding2core.port_b.Q_flow]
    @test sol[T_winding.T][end] >= 500 # not good but better than nothing
    @test sol[T_core.T] <= sol[T_winding.T]
end
