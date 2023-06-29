using ModelingToolkit, OrdinaryDiffEq #, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant

R = 1.0
C = 1.0
V = 1.0
@variables t
@named resistor = Resistor(R = R)
@named capacitor = Capacitor(C = C)
@named source = Voltage()
@named constant = Constant(k = V)
@named ground = Ground()

rc_eqs = [connect(constant.output, source.V)
    connect(source.p, resistor.p)
    connect(resistor.n, capacitor.p)
    connect(capacitor.n, source.n, ground.g)]

@named rc_model = ODESystem(rc_eqs, t,
    systems = [resistor, capacitor, constant, source, ground])
sys = structural_simplify(rc_model)
prob = ODEProblem(sys, [capacitor.v_start => 1], (0, 10.0))
sol = solve(prob, Tsit5())

plot(sol, vars = [capacitor.v, resistor.i],
    title = "RC Circuit Demonstration",
    labels = ["Capacitor Voltage" "Resistor Current"])
savefig("plot.png");

###### DC Motor

using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
# using Plots

@parameters t

R = 0.5 # [Ohm] armature resistance
L = 4.5e-3 # [H] armature inductance
k = 0.5 # [N.m/A] motor constant
J = 0.02 # [kg.m²] inertia
f = 0.01 # [N.m.s/rad] friction factor
tau_L_step = -0.3 # [N.m] amplitude of the load torque step

@named ground = Ground()
@named source = Voltage()
@named ref = Blocks.Step(height = 1, start_time = 0)
@named pi_controller = Blocks.LimPI(k = 1.1, T = 0.035, u_max = 10, Ta = 0.035)
@named feedback = Blocks.Feedback()
@named R1 = Resistor(R = R)
@named L1 = Inductor(L = L)
@named emf = EMF(k = k)
@named fixed = Fixed()
@named load = Torque(use_support = false)
@named load_step = Blocks.Step(height = tau_L_step, start_time = 3)
@named inertia = Inertia(J = J)
@named friction = Damper(d = f)
@named speed_sensor = SpeedSensor()

connections = [connect(fixed.flange, emf.support, friction.flange_b)
    connect(emf.flange, friction.flange_a, inertia.flange_a)
    connect(inertia.flange_b, load.flange)
    connect(inertia.flange_b, speed_sensor.flange)
    connect(load_step.output, load.tau)
    connect(ref.output, feedback.input1)
    connect(speed_sensor.w, :y, feedback.input2)
    connect(feedback.output, pi_controller.err_input)
    connect(pi_controller.ctr_output, :u, source.V)
    connect(source.p, R1.p)
    connect(R1.n, L1.p)
    connect(L1.n, emf.p)
    connect(emf.n, source.n, ground.g)]

@named model = ODESystem(connections, t,
    systems = [
        ground,
        ref,
        pi_controller,
        feedback,
        source,
        R1,
        L1,
        emf,
        fixed,
        load,
        load_step,
        inertia,
        friction,
        speed_sensor,
    ])
