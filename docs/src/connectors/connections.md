# Introduction
In Physical Network Acausal modeling each physical domain must define a **connector** to combine model components.  This is somewhat analogus to real life connections that are seen in electronics (i.e. battery connected to a wire) or fluid dynamics (i.e. pump connected to a pipe), to name a couple examples.  Each physical domain **connector** defines a minimum of 2 variables, one which is called a *Through* variable, and one which is called an *Across* variable.  Both Modelica and SimScape define these variables in the same way:

- [Modelica Connectors](https://mbe.modelica.university/components/connectors/#acausal-connection)
- [SimScape Connectors](https://www.mathworks.com/help/simscape/ug/basic-principles-of-modeling-physical-networks.html#bq89sba-6)

However, the standard libraries differ on the selection of the Across variable for the Mechanical Translation and Rotation libraries, Modelica choosing position and angle and SimScape choosing velocity and angular velocity, respectively for Translation and Rotation.  Modelica describes their decision [here](https://mbe.modelica.university/components/connectors/simple_domains/).  In summary they would like to provide less integration in the model to avoid lossy numerical behavior, but this decision assumes the lowest order derivative is needed by the model.  Numerically it is possible to define the connector either way, but there are some consequences to this decision, and therefore we will study them in detail here as they relate to ModelingToolkit.  

# Through and Across Variable Theory
### General
The idea behind the selection of the **through** variable is that it should be a time derivative of some conserved quantity. The conserved quantity should be expressed by the **across** variable.  In general terms the physical system is given by

- Energy Dissipation & Flow: 

```math
\begin{aligned}
    \partial {\color{blue}{across}} / \partial t \cdot c_1 = {\color{green}{through}}  \\
    {\color{green}{through}} \cdot c_2 = {\color{blue}{across}} 
\end{aligned}
```



### Electrical
So for the Electrical domain the across variable is *voltage* and the through variable *current*.  Therefore 

- Energy Dissipation: 
```math
\partial {\color{blue}{voltage}} / \partial t \cdot capacitance = {\color{green}{current}} 
```

- Flow: 
```math
{\color{green}{current}} \cdot resistance = {\color{blue}{voltage}}
```

### Translational
And for the translation domain, choosing *velocity* for the across variable and *force* for the through gives

- Energy Dissipation: 
```math  
\partial {\color{blue}{velocity}} / \partial t \cdot mass = {\color{green}{force}} 
```

- Flow: 
```math 
{\color{green}{force}} \cdot (1/damping) = {\color{blue}{velocity}} 
```

The diagram here shows the similarity of problems in different physical domains.  

![Through and Across Variables](through_across.png)

### Translational using Position Across Variable
Now, if we choose *position* for the across variable, a similar relationship can be established, but the patern must be broken.

- Energy Dissipation: 
```math  
\partial^2 {\color{blue}{position}} / \partial t^2 \cdot mass = {\color{green}{force}} 
```

- Flow: 
```math 
{\color{green}{force}} \cdot (1/damping) = \partial {\color{blue}{position}} / \partial t 
```

As can be seen, we must now establish a higher order derivative to define the Energy Dissipation and Flow equations, requiring an extra equation, as will be shown in the example below.

# Examples
### Electrical Domain
We can generate the above relationship with ModelingToolkit and the ModelingToolkitStandardLibrary using 3 blocks:

- Capacitor: for energy storage with initial voltage = 1V
- Resistor: for energy flow
- Ground: for energy sink

As can be seen, this will give a 1 equation model matching our energy dissipation relationship

```@example connections
using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq
using Plots

@parameters t

@named resistor = Resistor(R = 1)
@named capacitor = Capacitor(C = 1)
@named ground = Ground()

eqs = [
    connect(capacitor.p, resistor.p)
    connect(resistor.n, ground.g, capacitor.n)
    ]

@named model = ODESystem(eqs, t; systems=[resistor, capacitor, ground])

sys = structural_simplify(model)

println.(equations(sys))
nothing # hide
```

The solution shows what we would expect, a non-linear disipation of voltage and releated decrease in current flow...

```@example connections
prob = ODEProblem(sys, [1.0], (0, 10.0), [])
sol = solve(prob, ImplicitMidpoint(); dt=0.01)

p1=plot(sol, idxs=[capacitor.v])
p2=plot(sol, idxs=[resistor.i])
plot(p1,p2)
savefig("electrical.png"); nothing # hide
```

![Plot of Electrical Example](electrical.png)

### Mechanical Translational Domain 
#### Across Variable = velocity
Now using the Translational library based on velocity, we can see the same relationship with a system reduced to a single equation, using the components:

- Body (i.e. moving mass): for kinetic energy storage with an initial velocity = 1m/s
- Damper: for energy flow
- Fixed: for energy sink

```@example connections
using ModelingToolkitStandardLibrary
const TV = ModelingToolkitStandardLibrary.Mechanical.Translational

@named damping = TV.Damper(d=1, v1_0=1)
@named body = TV.Mass(m=1, v_0=1)
@named ground = TV.Fixed()

eqs = [
    connect(damping.port_a, body.port)
    connect(ground.port, damping.port_b)
    ]

@named model = ODESystem(eqs, t; systems=[damping, body, ground])

sys = structural_simplify(model)

println.(full_equations(sys))
nothing # hide
```

As expected we have a similar solution...
```@example connections
prob = ODEProblem(sys, [], (0, 10.0), [])
sol_v = solve(prob, ImplicitMidpoint(); dt=0.01)

p1=plot(sol_v, idxs=[body.v])
p2=plot(sol_v, idxs=[damping.f])
plot(p1,p2)
savefig("mechanical_velocity.png"); nothing # hide
```

![Plot of Mechanical (Velocity Based) Example](mechanical_velocity.png)


#### Across Variable = position
Now, let's consider the position based approach.  We can build the same model with the same components.  As can be seen, we now end of up with 2 equations, because we need to relate the lower derivative (position) to force (with acceleration).  

```@example connections
const TP = ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition

@named damping = TP.Damper(d=1, v1_0=1)
@named body =    TP.Mass(m=1, v_0=1)
@named ground =  TP.Fixed(s_0=0)

    eqs = [
        connect(damping.port_a, body.port)
        connect(ground.port, damping.port_b)
        ]

@named model = ODESystem(eqs, t; systems=[damping, body, ground])

sys = structural_simplify(model)

println.(full_equations(sys))
nothing # hide
```
As can be seen, we get exactly the same result.  The only difference here is that we are solving an extra equation, which allows us to plot the body position as well.

```@example connections
prob = ODEProblem(sys, [], (0, 10.0), [])
sol_p = solve(prob, ImplicitMidpoint(); dt=0.01)

p1=plot(sol_p, idxs=[body.v])
p2=plot(sol_p, idxs=[damping.f])
p3=plot(sol_p, idxs=[body.s])

plot(p1,p2,p3)
savefig("mechanical_position.png"); nothing # hide
```

![Plot of Mechanical (Velocity Based) Example](mechanical_position.png)

The question then arises, can the position be plotted when using the Mechanical Translational Domain based on the Velocity Across variable?  Yes, we can!  There are 2 solutions:

1. the `Mass` component will add the position variable when the `s_0` parameter is used to set an initial position.  Otherwise the position is not tracked by the component.

```julia
@named body = TV.Mass(m=1,  v_0=1,  s_0=0)
```

2. implement a `PositionSensor`
TODO: Implement Translation Sensors


Either option will produce the same result regardless of which across variable is used.  If the same result is given, why are both options included in the Standard Library, what are the differences?  These differences will be discussed next so that an informed decision can be made about which domain is best for your model.

# Mechanical/Translational Library Differences (Velocity vs. Position Connectors)
## Initialization
The main difference between `ModelingToolkitStandardLibrary.Mechanical.Translational` and `ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition` is how they are initialized.  In the `ModelingToolkitStandardLibrary` initialization parameters are defined at the component level, so we simply need to be careful to set the correct initial conditions for the domain that it used.  Let's use the following example problem to explain the differences.  

![Example Mechanical Model](model.png)

In this problem we have a mass, spring, and damper which are connected to a fixed point.  Let's see how each component is defined.

#### Damper
The damper will connect the flange/port 1 (`port_a`) to the mass, and flange/port 2 (`port_b`) to the fixed point.  For both position and velocity based domains, we set the damping constant `d=1` and `v1_0=1` and leave the default for `v2_0` at 0.  For the position domain we also need to set the initial positions for `port_a` and `port_b`.

```@example connections
@named dv = TV.Damper(d=1, v1_0=1)
@named dp = TP.Damper(d=1, v1_0=1, s1_0=3, s2_0=1)
nothing # hide
```

#### Spring
The spring will connect the flange/port 1 (`port_a`) to the mass, and flange/port 2 (`port_b`) to the fixed point.  For both position and velocity based domains, we set the spring constant `k=1`.  The velocity domain then requires the initial velocity `v1_0` and initial spring stretch `delta_s_0`.  The position domain instead needs the initial positions for `port_a` and `port_b` and the natural spring length `l`.  

```@example connections
@named sv = TV.Spring(k=1, v1_0=1, delta_s_0=1)
@named sp = TP.Spring(k=1, s1_0=3, s2_0=1, l=1)
nothing # hide
```

#### Mass
For both position and velocity based domains, we set the mass `m=1` and initial velocity `v_0=1`. Like the damper, the position domain requires the position initial conditions set as well.  

```@example connections
@named bv = TV.Mass(m=1, v_0=1)
@named bp = TP.Mass(m=1, v_0=1, s_0=3)
nothing # hide
```

#### Fixed
Here the velocity domain requires no initial condition, but for our model to work as defined we must set the position domain component to the correct intital position.

```@example connections
@named gv =  TV.Fixed()
@named gp =  TP.Fixed(s_0=1)
nothing # hide
```

### Comparison
As can be seen, the position based domain requires more initial condition information to be properly defined since the absolute position information is required.  Thereore based on the model being described, it may be more natural to choose one domain over the other.  

Let's define a quick function to simplify and solve the 2 different systems.

```@example connections
function simplify_and_solve(damping, spring, body, ground)

       eqs = [connect(spring.port_a, body.port, damping.port_a)
              connect(spring.port_b, damping.port_b, ground.port)
              ]

       @named model = ODESystem(eqs, t; systems = [ground, body, spring, damping])

       sys = structural_simplify(model)

       println.(full_equations(sys))

       prob = ODEProblem(sys, [], (0, 10.0), [])
       sol = solve(prob, ImplicitMidpoint(), dt=0.01)

       return sol,sys
end
nothing # hide
```

Now let's solve the velocity domain model

```@example connections
solv,sysv=simplify_and_solve(dv, sv, bv, gv);
nothing # hide
```

And the position domain model

```@example connections
solp,sysp=simplify_and_solve(dp, sp, bp, gp);
nothing # hide
```

Now we can plot the comparison of the 2 models and see they give the same result.  
```@example connections
plot(ylabel="mass velocity [m/s]")
plot!(solv, idxs=[bv.v])
plot!(solp, idxs=[bp.v])
savefig("mass_velocity.png"); nothing # hide
```

![Mass Velocity Comparison](mass_velocity.png)


But, what if we wanted to plot the mass position?  This is easy for the position based domain, we have the state `bp₊s(t)`, but for the velocity based domain we have `sv₊delta_s(t)` which is the spring stretch.  To get the absolute position we add the spring natrual length (1m) and the fixed position (1m).  As can be seen, we then get the same result.

```@example connections
plot(ylabel="mass position [m]")
plot!(solv.t, solv[sv.delta_s] .+ 1 .+ 1, label="sv.delta_s(t) + 1 + 1")
plot!(solp, idxs=[bp.s])
savefig("mass_position.png"); nothing # hide
```

![Mass Position Comparison](mass_position.png)

So in conclusion, the position based domain gives easier access to absolute position information, but requires more initial condition information.  


## Acuracy
One may ask then what is the trade off in terms of numerical acuracy?  When we look at the simplified equations, we can see that actually both systems solve the same equations.  The differential equations of the velocity domain are

```math
\begin{aligned}
m \cdot \dot{v} +  d \cdot v + k \cdot delta_s = 0  \\
\dot{delta_s} = v
\end{aligned}
```

And for the position domain are

```math
\begin{aligned}
m \cdot \dot{v} +  d \cdot v + k \cdot (s - s2_0 - l) = 0   \\
\dot{s} = v
\end{aligned}
```

By definition the spring stretch is

```math
delta_s = s - s2_0 - l
```

Which means both systems are actually solving the same exact system.  We can plot the numerical difference between the 2 systems and see the result is negligible.

```@example connections
plot(title="numerical difference: vel. vs. pos. domain", xlabel="time [s]", ylabel="solv[bv.v] .- solp[bp.v]")
plot!(solv.t, solv[bv.v] .- solp[bp.v], label="")
savefig("err.png"); nothing # hide
```

![Accuracy Comparison](err.png)

