using ModelingToolkitStandardLibrary.Mechanical.Translational, ModelingToolkit, OrdinaryDiffEq,
      Test
import ModelingToolkitStandardLibrary.Blocks


@parameters t
D = Differential(t)

@testset "one falling body" begin
    #body starts at 1m, moving up at 2m/s in gravity of 10
    @named body = Body(x0 = 1.0, v0 = 2.0, g=-10.0)

    @named model = ODESystem(body.T.f ~ 0, t, [], []; systems=[body])

    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0, 1.0))
    sol = solve(prob, ImplicitMidpoint(), dt=0.1)

    #x = g*t^2/2 + v_int*t + x_int
    #dx = g*t + v_int

    @test sol.retcode == :Success
    @test sol[body.x, 1] == 1.0
    @test sol[body.dx, 1] == 2.0
    @test sol[body.x, end] ≈ -10/2 + 2.0 + 1.0

#=
    using CairoMakie
    lines(sol.t, sol[body.x])
=#

end


@testset "two bodies ForcedBoundary" begin
    @named body1 = Body(ForcedBoundary; mass=100.0, x0 = 10.0, g=-10.0, lower_limit = 0.0)
    @named body2 = Body(mass = 1000.0, g = -10.0)
    @named spring = Spring(k=1000.0)

    eqs = [
        connect(body1.T, spring.T1)
        connect(body2.T, spring.T2)
        ]

    @named model = ODESystem(eqs, t, [], []; systems=[body1, body2, spring])

    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0, 20.0))
    
    #
    
    sol = solve(prob, ImplicitEuler(nlsolve = NLNewton(check_div = false, always_new = true, max_iter=100, relax=0.0)),  
                dt = 4e-4,
                adaptive=false)

    @test sol(5.0; idxs=body1.x) ≈ 0.0 atol=1e-7
end

#=
    using CairoMakie
    begin
        f = Figure()
        a = Axis(f[1,1])
        lines!(a, sol.t, sol[body1.x])   
        lines!(a, sol.t, sol[body2.x])
        # ylims!(a, -0.0001, 0.0001)

        
        a = Axis(f[2,1])
        lines!(a, sol.t, sol[body1.f_ceil])
        lines!(a, sol.t, sol[body1.f_floor])
        lines!(a, sol.t, sol[spring.f])
        # ylims!(a, -100000, 0)
        
        a = Axis(f[3,1])
        lines!(a, sol.t, sol[body1.x] .== 0)
        
        f   
    end

    using Makie
    tm = Observable(0.0)
    begin
        f = Figure()
        a = Axis(f[1,2], aspect=DataAspect(), )
        hidedecorations!(a)
        s = @lift(sol($tm))

        m1x = @lift($s[2])
        m2x = @lift($s[4])

        sz1 = 0.5
        lines!(a, [-sz1, sz1, sz1, -sz1, -sz1], @lift([$m1x, $m1x, $m1x+sz1*2, $m1x+sz1*2, $m1x]))

        sz = 1.0
        lines!(a, [-sz, sz, sz, -sz, -sz], @lift([$m2x, $m2x, $m2x+sz*2, $m2x+sz*2, $m2x]))

        lines!(a, [0.0, 0.0], @lift([$m1x+0.5, $m2x+1.0]), color=:red, width=2, linestyle=:dot)

        hlines!(a, 0.0, color=:gray)
        hlines!(a, -10.0, color=:gray, linestyle=:dash)

        ylims!(a, -40, 15)

        a = Axis(f[1, 1], xlabel="time [s]", ylabel="position [m]")
        lines!(a, sol.t, sol[body1.x])   
        lines!(a, sol.t, sol[body2.x])

        scatter!(a, tm, m1x)
        scatter!(a, tm, m2x)
        ylims!(a, -40, 15)

        f
    end

    framerate = 30
    timestamps = range(0, 20, step=1/framerate)

    record(f, "time_animation.mp4", timestamps;
            framerate = framerate) do t
        tm[] = t
    end
=#   



