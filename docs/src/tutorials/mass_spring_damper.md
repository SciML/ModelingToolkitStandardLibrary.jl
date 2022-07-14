# Translational Mass Spring Damper System

This tutorial demonstrates the response of a mass-spring-damper system. Rather than applying an external force, the mass is initialized with some positional offset from equilibrium and some non-zero velocity (system responses under such conditions are commonly referred to as "unforced" responses). Assuming a lossless system, a mass-spring system in the absence of a damper will continually oscillate at some natural frequency. The addition of a damper will result in (1) underdamped, (2) overdamped, or (3) critically damped motion, depending on the relative values of our mass (`m`), spring constant (`c`), and damping constant (`d`). In this first example, we will demonstrate an underdamped system. Values for `m`, `c`, and `d` were chosen such that the roots of the polynomial `s^2 + d/m * s + c/m` are complex conjugates (the equivalent condition for underdamped motion is `d < sqrt(4 * m * c)`). In the case of underdamped systems, our mass will undergo damped oscillations. This damping corresponds to system energy loss, and the amplitudes of said oscillations will decrease towards zero as time approaches infinity.

Let's start by importing what we need.

```@example
using ModelingToolkitStandardLibrary.Mechanical.Translational, ModelingToolkit, OrdinaryDiffEq, ControlSystems, Plots
import ModelingToolkitStandardLibrary.Blocks
```

Now let's initialize our system state and define our components and connections. Here, `s0` is the fixed offset position of our flange housing, `s_rel0` is our unstretched spring length, and `s_start` and `v_start` are the initial values of the absolute position and linear velocity of our sliding mass, respectively. Finally, `m`, `c`, and `d` are our mass, spring, and damping constants, respectively. All units are in SI units.

```@parameters t
D = Differential(t)

s0 = 4.5
s_rel0 = 1
s_start = 3
v_start = 10
m = 1
c = 10
d = 1

@named fixed = Fixed(s0=s0)
@named mass = Mass(m=m, s_start=s_start, v_start=v_start)
@named damper = Damper(d=d)
@named spring = Spring(c=c, s_rel0=s_rel0)

connections = [
    connect(mass.flange_b, damper.flange_a, spring.flange_a)
    connect(damper.flange_b, spring.flange_b, fixed.flange)
]
```

Before solving our ODE system, let us ensure that the condition for underdamped motion is met.

```
@assert(d < sqrt(4 * m * c))
```

Now let's solve our ODE system and plot our results (note that we are plotting the position of our mass as a function of time).

```@named model = ODESystem(connections, t, systems=[fixed, mass, damper, spring])
sys = structural_simplify(model)
prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 10.0), saveat=0.002)
sol = solve(prob, DFBDF(), abstol=1e-12, reltol=1e-12)

plot(sol, vars = [mass.s], 
    title = "Mass Spring Damper - Underdamped Motion", 
    xlabel = "Time (s)", 
    ylabel = "Position (m)")

savefig("mass_spring_damper_underdamped.png"); nothing # hide
```

!["Position vs. Time Plot of a Mass in an Underdamped Mass-Spring-Damper System](mass_spring_damper_underdamped.png)

We'll now demonstrate a mass-spring-damper system that is overdamped. For such systems, `d > sqrt(4 * m * c)`, and the zeros of our characteristic polynomial are real. Let us re-initialize our `m`, `c`, and `d` constants accordingly, and re-define our components and connections.

```
m, c, d = 1, 1, 10

@named mass = Mass(m=m, s_start=s_start, v_start=v_start)
@named damper = Damper(d=d)
@named spring = Spring(c=c, s_rel0=s_rel0)

connections = [
    connect(mass.flange_b, damper.flange_a, spring.flange_a)
    connect(damper.flange_b, spring.flange_b, fixed.flange)
]
```

Let us ensure that the condition for overdamped motion is met.

```
@assert(d > sqrt(4 * m * c))
```

Let's solve our ODE system and plot our results. Observe that for overdamped systems, motion is non-oscillatory.

```@named model = ODESystem(connections, t, systems=[fixed, mass, damper, spring])
sys = structural_simplify(model)
prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 10.0), saveat=0.002)
sol = solve(prob, DFBDF(), abstol=1e-12, reltol=1e-12)

plot(sol, vars = [mass.s], 
    title = "Mass Spring Damper - Overdamped Motion", 
    xlabel = "Time (s)", 
    ylabel = "Position (m)")

savefig("mass_spring_damper_overdamped.png"); nothing # hide
```

!["Position vs. Time Plot of a Mass in an Overdamped Mass-Spring-Damper System](mass_spring_damper_overdamped.png)

Our last demonstration will be of a mass-spring-damper system that is critically damped. The condition for critically damped motion is `d = sqrt(4 * m * c)`, where the roots of our characteristic polynomial are both equal to `d/2m`. Again, let us re-initialize and re-define our system state, components, and connections, solve our ODE system, and plot our results. Observe that for critically damped systems, motion is non-oscillatory.

```
m, c = 1, 1
d = sqrt(4 * m * c)

@named mass = Mass(m=m, s_start=s_start, v_start=v_start)
@named damper = Damper(d=d)
@named spring = Spring(c=c, s_rel0=s_rel0)

connections = [
    connect(mass.flange_b, damper.flange_a, spring.flange_a)
    connect(damper.flange_b, spring.flange_b, fixed.flange)
]

@named model = ODESystem(connections, t, systems=[fixed, mass, damper, spring])
sys = structural_simplify(model)
prob = DAEProblem(sys, D.(states(sys)) .=> 0.0, [D(D(mass.s)) => 1.0], (0, 10.0), saveat=0.002)
sol = solve(prob, DFBDF(), abstol=1e-12, reltol=1e-12)

plot(sol, vars = [mass.s], 
    title = "Mass Spring Damper - Critically Damped Motion", 
    xlabel = "Time (s)", 
    ylabel = "Position (m)")

savefig("mass_spring_damper_critically_damped.png"); nothing # hide
```

!["Position vs. Time Plot of a Mass in a Critically Damped Mass-Spring-Damper System](mass_spring_damper_critically_damped.png)