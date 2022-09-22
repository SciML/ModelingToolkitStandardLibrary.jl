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
    \partial {\color{blue}{across}} / \partial t \cdot c_1 = {\color{green}{through}} 
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
\color{green}{current} \cdot resistance = \color{blue}{voltage}
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
    connect(capacitor.p, resistor.n)
    connect(resistor.p, ground.g, capacitor.n)
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
module TranslationalVelocity
    using ModelingToolkit
    using ModelingToolkitStandardLibrary.Mechanical.Translational

    @parameters t

    @named damping = Damper(d = 1)
    @named body = Body(m = 1, v0=1)
    @named ground = Fixed()

    eqs = [
        connect(damping.port_a, body.port)
        connect(ground.port, damping.port_b)
        ]

    @named model = ODESystem(eqs, t; systems=[damping, body, ground])

    sys = structural_simplify(model)
end

sys = TranslationalVelocity.sys
println.(full_equations(sys))
nothing # hide
```

As expected we have a similar solution...
```@example connections
prob = ODEProblem(sys, [1.0], (0, 10.0), [])
sol_v = solve(prob, ImplicitMidpoint(); dt=0.01)

p1=plot(sol_v, idxs=[TranslationalVelocity.body.v])
p2=plot(sol_v, idxs=[TranslationalVelocity.damping.f])
plot(p1,p2)
savefig("mechanical_velocity.png"); nothing # hide
```

![Plot of Mechanical (Velocity Based) Example](mechanical_velocity.png)


#### Across Variable = position
Now, let's consider the position based approach.  We can build the same model with the same components.  As can be seen, we now end of up with 2 equations, because we need to relate the lower derivative (position) to force (with acceleration).  

```@example connections
module TranslationalPosition
    using ModelingToolkit
    using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition

    @parameters t
    D = Differential(t)

    # Let's define a simple body that only tracks the across and through variables...
    function Body(; name, m, v0 = 0.0)
        @named flange = Flange()
        pars = @parameters m=m v0=v0
        vars = @variables begin
            v(t) = v0
            f(t) = m*v0
        end
        eqs = [
            D(flange.s) ~ v
            flange.f ~ f

            D(v) ~ f/m
            ]
        return compose(ODESystem(eqs, t, vars, pars; name = name), flange)
    end

    @named damping = Damper(d = 1)
    @named body = Body(m = 1, v0=1)
    @named ground = Fixed()

    eqs = [
        connect(damping.flange_b, body.flange)
        connect(ground.flange, damping.flange_a)
        ]

    @named model = ODESystem(eqs, t; systems=[damping, body, ground])

    sys = structural_simplify(model)
end

sys = TranslationalPosition.sys

println.(full_equations(sys))
nothing # hide
```

As can be seen, we get exactly the same result.  The only difference here is that we are solving an extra equation, which allows us to plot the body position as well.

```@example connections
prob = ODEProblem(sys, [0.0, 1.0], (0, 10.0), [])
sol_p = solve(prob, ImplicitMidpoint(); dt=0.01)

p1=plot(sol_p, idxs=[TranslationalPosition.body.v])
p2=plot(sol_p, idxs=[TranslationalPosition.damping.f])
p3=plot(sol_p, idxs=[TranslationalPosition.body.flange.s])

plot(p1,p2,p3)
savefig("mechanical_position.png"); nothing # hide
```

![Plot of Mechanical (Velocity Based) Example](mechanical_position.png)

The question then arises, can the position be plotted when using the Mechanical Translational Domain based on the Velocity Across variable?  Yes, we can!  There are 2 solutions:

1. the `Body` component will add the position variable when the `s0` parameter is used to set an initial position.  Otherwise the position is not tracked by the component.

```julia
@named body = Body(m = 1, v0=1, s0=0)
```

2. implement a `PositionSensor`

```julia
@named damping = Damper(d = 1)
@named body = Body(m = 1, v0=1)
@named ground = Fixed()
@named sensor = PositionSensor(s0=0)

eqs = [
    connect(damping.port_a, body.port, sensor.port)
    connect(ground.port, damping.port_b)
    ]
```

Either option will produce the same result as shown for the Mechanical Translational Domain based on the Position Across variable.  If the same result is given, why are both options included in the Standard Library, what are the differences?  These differences will be discussed next so that an informed decision can be made about which domain is best for your model.

# Differences
## Initialization

```julia
module TranslationalVelocity
    using ModelingToolkit
    using ModelingToolkitStandardLibrary.Mechanical.Translational

    @parameters t

    @named damping = Damper(d=1)
    @named spring = Spring(k=1)
    @named body = Body(m = 1, v0=1)
    @named ground = Fixed()

    eqs = [
        connect(damping.port_a, spring.port_a, body.port)
        connect(ground.port, damping.port_b, spring.port_b)
        ]

    @named model = ODESystem(eqs, t; systems=[damping, body, ground, spring])

    sys = structural_simplify(model)
end
sys = TranslationalVelocity.sys
prob = ODEProblem(sys, [1.0, 0.0], (0, 10.0), [])
sol = solve(prob, ImplicitMidpoint(); dt=0.01)
```

## Acuracy





