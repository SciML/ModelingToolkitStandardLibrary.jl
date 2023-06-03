using ModelingToolkit
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.Translational
using DifferentialEquations
using Setfield
using Test

@parameters t

@named link1 = Link(; m = 1, l = 10, I = 85, g = -9.807)
@named link2 = Link(; m = 1, l = 10, I = 85, g = -9.807, x1_0=10)
@named link3 = Link(; m = 1, l = 10, I = 85, g = -9.807, x1_0=20)
@named link4 = Link(; m = 1, l = 10, I = 85, g = -9.807, x1_0=30)
@named link5 = Link(; m = 1, l = 10, I = 85, g = -9.807, x1_0=40)

@named joint1 = RevoluteJoint(d=1)
@named joint2 = RevoluteJoint(d=200)
@named joint3 = RevoluteJoint(d=200)
@named joint4 = RevoluteJoint(d=200)
@named joint5 = RevoluteJoint(d=200)

@named cart = Mass(; m = 5, s_0 = 0)
@named mb2t = MultiBody2Translational()

eqs = [

    connect(cart.flange, mb2t.T)

    connect(mb2t.M,   joint1.M1)
    connect(link1.M1, joint1.M2)
    connect(link1.M2, joint2.M1)
    connect(link2.M1, joint2.M2)
    connect(link2.M2, joint3.M1)
    connect(link3.M1, joint3.M2)
    connect(link3.M2, joint4.M1)
    connect(link4.M1, joint4.M2)
    connect(link4.M2, joint5.M1)
    connect(link5.M1, joint5.M2)
    # connect_mbs(link5.M2)

       ]

@named model = ODESystem(eqs, t, [], []; systems = [cart, mb2t, link1, link2, link3, link4, link5, joint1, joint2, joint3, joint4, joint5]) 

sys = structural_simplify(model)
unset_vars = setdiff(states(sys), keys(ModelingToolkit.defaults(sys)))
prob = ODEProblem(sys, unset_vars .=> 0.0, (0.0, 20), []; jac = true)
NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 4 // 10)
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON), dt = 4e-4, adaptive = false)

@test sol[link1.x1][end] ≈ 18.48185172484587

#=
using CairoMakie
begin
    f = Figure()
    a = Axis(f[1,1],xlabel="time [s]", ylabel="cart x pos. [m]")
    lines!(a, sol.t, sol[link1.x1])
    f
end

function plot_link(sol, sys, tmax)
    tm = Observable(0.0)
    idx = Dict(reverse.(enumerate(states(sys))))

    fig = Figure()
    a = Axis(fig[1,1], aspect=DataAspect(), )
    hidedecorations!(a)
    s = @lift(sol($tm, idxs=[link1.x1, link1.x2, link1.y1, link1.y2, link2.x1, link2.x2, link2.y1, link2.y2, link3.x1, link3.x2, link3.y1, link3.y2, link4.x1, link4.x2, link4.y1, link4.y2, link5.x1, link5.x2, link5.y1, link5.y2]))

    m1x1 = @lift($s[1])
    m1x2 = @lift($s[2])
    m1y1 = @lift($s[3])
    m1y2 = @lift($s[4])

    m2x1 = @lift($s[5])
    m2x2 = @lift($s[6])
    m2y1 = @lift($s[7])
    m2y2 = @lift($s[8])

    m3x1 = @lift($s[9])
    m3x2 = @lift($s[10])
    m3y1 = @lift($s[11])
    m3y2 = @lift($s[12])

    m4x1 = @lift($s[13])
    m4x2 = @lift($s[14])
    m4y1 = @lift($s[15])
    m4y2 = @lift($s[16])

    m5x1 = @lift($s[17])
    m5x2 = @lift($s[18])
    m5y1 = @lift($s[19])
    m5y2 = @lift($s[20])

    sz1 = 0.5
    # lines!(a, [-sz1, sz1, sz1, -sz1, -sz1], @lift([$m1x1, $m1x1, $m1x1+sz1*2, $m1x1+sz1*2, $m1x1]))

    hlines!(a, 0, color=:gray)
    # scatter!(a, m1x1[], marker=:square, markersize=10, color=:black)
    lines!(a, @lift([$m1x1, $m1x2]) , @lift([$m1y1, $m1y2]), linewidth=10, color=:blue)
    lines!(a, @lift([$m2x1, $m2x2]) , @lift([$m2y1, $m2y2]), linewidth=10, color=:red)
    lines!(a, @lift([$m3x1, $m3x2]) , @lift([$m3y1, $m3y2]), linewidth=10, color=:blue)
    lines!(a, @lift([$m4x1, $m4x2]) , @lift([$m4y1, $m4y2]), linewidth=10, color=:red)
    lines!(a, @lift([$m5x1, $m5x2]) , @lift([$m5y1, $m5y2]), linewidth=10, color=:blue)

    CairoMakie.ylims!(a, -50, 20)
    CairoMakie.xlims!(a, -20, 50)

    # a = Axis(fig[1, 1], xlabel="time [s]", ylabel="position [m]")
    # lines!(a, sol.t, sol[2,:])
    # lines!(a, sol.t, sol[4,:])

    # scatter!(a, tm, m1x)
    # scatter!(a, tm, m2x)
    # ylims!(a, -60, 30)

    framerate = 30
    timestamps = range(0, tmax, step=1/framerate)

    record(fig, "links.mp4", timestamps;
            framerate = framerate) do t
        tm[] = t
    end

    #=
    CairoMakie.Makie.Record(fig, timestamps; framerate=framerate) do t
        tm[] = t
    end
    =#

    nothing
end

plot_link(sol, sys, 20)

=#