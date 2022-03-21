using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks: t, Saturation, DeadZone, Integrator

#=
Testing strategy:
The general strategy is to test systems using simple intputs where the solution is known on closed form. For algebraic systems (without differential variables), an integrator with a constant input is often used together with the system under test. 
=#

## Saturation
@testset "Saturation" begin
    @info "Testing Saturation"
    y_max = 0.8
    y_min = -0.6
    @named sat = Saturation(; y_min, y_max)
    @test count(ModelingToolkit.isinput, states(sat)) == 1
    @test count(ModelingToolkit.isoutput, states(sat)) == 1
    @named iosys = ODESystem([sat.u~sin(t)], t, systems=[sat])

    sys = structural_simplify(iosys)
    stateind(sym) = findfirst(isequal(sym),states(sys))
    prob = ODEProblem(sys, Pair[sat.u=>0], (0.0, 2pi))
    sol = solve(prob, Rosenbrock23(), saveat=0:0.1:2pi)
    y = clamp.(sin.(sol.t), y_min, y_max)
    # plot([sol[sat.y] y])
    @test sol[sat.y] ≈ y rtol = 1e-6
end


@testset "DeadZone" begin
    @info "Testing DeadZone"

    ie = ifelse
    deadzone(u, u_min, u_max) = ie(u > u_max, u-u_max, ie( u < u_min, u-u_min, 0))

    u_max = 1
    u_min = -2
    @named dz = DeadZone(; u_min, u_max)
    @named int = Integrator()
    @test count(ModelingToolkit.isinput, states(dz)) == 1
    @test count(ModelingToolkit.isoutput, states(dz)) == 1
    @named iosys = ODESystem([
        int.u ~ 1
        int.y ~ dz.u
        ], t, systems=[dz, int]
    )

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => -3], (0.0, 5))
    sol = solve(prob, Rosenbrock23(), saveat=0:0.1:5)
    y = deadzone.(sol[int.y], u_min, u_max)
    @test y[1] == -3 - u_min
    @test y[end] ≈ 2 - u_max
    @test y[end÷2] == 0
    
    # plot([y sol[dz.y] sol[dz.u]])
    @test sol[dz.y] ≈ y rtol = 1e-6
end