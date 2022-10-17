# DC Motor with PI-controller
In this example a PI-controller is setup for speed control of a DC-motor. An equivalent circuit diagram is depicted below.

![DC-motor](https://user-images.githubusercontent.com/50108075/196108356-0e8605e3-61a9-4006-8559-786252e55928.png)

The electrical part consists of a resistance and inductance. The coupling between the electrical and rotational domain is done via an electro motive force (EMF) component. The voltage across the EMF is proportional to the angular velocity and the current is proportional to the torque. On the mechanical side viscous friction in e.g. a bearing and the inertia of the shaft is modelled.

A PI-controller with anti windup measure should be used as a speed controller. A simulation is performed to verify the tracking performance of the controller and the disturbance rejection capabilities. 

First the needed packages are imported and the parameters of the model defined.

```@example dc_motor_pi
using ModelingToolkit
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using Plots

@parameters t

R = 0.5 # [Ohm] armature resistance
L = 4.5e-3 # [H] armature inductance
k = 0.5 # [N.m/A] motor constant
J = 0.02 # [kg.m²] inertia
f = 0.01 # [N.m.s/rad] friction factor
tau_L_step = -0.3 # [N.m] amplitude of the load torque step 
nothing # hide
```

The actual model can now be composed.

```@example dc_motor_pi
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
               connect(speed_sensor.w, feedback.input2)
               connect(feedback.output, pi_controller.err_input)
               connect(pi_controller.ctr_output, source.V)
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
nothing # hide
```

Now the model can be simulated. Typical rotational mechanical systems are described via `DAE` 
(differential algebraic equations), however in this case ModelingToolkit can simplify the model enough
so that it can be represented as a system of `ODEs` (ordinary differential equations).

```@example dc_motor_pi
sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 6.0))
sol = solve(prob, Rodas4())

p1 = Plots.plot(sol.t, sol[inertia.w], ylabel = "Angular Vel. in rad/s",
                label = "Measurement", title = "DC Motor with Speed Controller")
Plots.plot!(sol.t, sol[ref.output.u], label = "Reference")
p2 = Plots.plot(sol.t, sol[load.tau.u], ylabel = "Disturbance in Nm", label = "")
Plots.plot(p1, p2, layout = (2, 1))
```
