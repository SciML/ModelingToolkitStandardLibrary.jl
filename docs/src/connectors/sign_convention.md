# Sign Convention

```@example sign_convention
using ModelingToolkit

@parameters t
D = Differential(t)

@mtkmodel System begin
    @components begin
        mass = Mass(; m=10)
        force = ConstantForce2()
    end
    @equations begin
        connect(mass.flange, force.flange)
    end
end
@mtkbuild sys = System()
foreach(println, full_equations(sys))
```

```@example sign_convention
using ModelingToolkitStandardLibrary.Electrical

@mtkmodel System begin
    @components begin
        capacitor = Capacitor(; C=10)
        current = ConstantCurrent()
        ground = Ground()
    end
    @equations begin
        connect(current.p, capacitor.n)
        connect(capacitor.p, ground.g, current.n)
    end
end
@mtkbuild sys = System()
foreach(println, full_equations(sys))


@mtkmodel System begin
    @components begin
        capacitor = Capacitor(; C=10)
        current = ConstantCurrent()
        ground = Ground()
    end
    @equations begin
        connect(current.n, capacitor.p)
        connect(capacitor.n, ground.g, current.p)
    end
end
@mtkbuild sys = System()
foreach(println, full_equations(sys))
```


```@example sign_convention
using ModelingToolkitStandardLibrary.Electrical

@mtkmodel System begin
    @components begin
        capacitor = Capacitor(; C=10)
        current = ConstantCurrent()
        ground = Ground()
    end
    @equations begin
        connect(current.p, capacitor.n)
        connect(capacitor.p, ground.g, current.n)
    end
end
@mtkbuild sys = System()
foreach(println, full_equations(sys))


@mtkmodel System begin
    @components begin
        capacitor = Capacitor(; C=10)
        current = ConstantCurrent()
        ground = Ground()
    end
    @equations begin
        connect(current.n, capacitor.p)
        connect(capacitor.n, ground.g, current.p)
    end
end
@mtkbuild sys = System()
foreach(println, full_equations(sys))
```


```@example sign_convention
using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible

@mtkmodel System begin
    @components begin
        volume = FixedVolume(; vol=10.0, p_int=0.0)
        flow = ConstantMassFlow()
        fluid = HydraulicFluid()
    end
    @equations begin
        connect(flow.port, volume.port)
        connect(fluid, flow.port)
    end
end
@mtkbuild sys = System()
foreach(println, full_equations(sys))
```