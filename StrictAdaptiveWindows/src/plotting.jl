import Plots

import .._!
import ..calculate_timegrid
import ..calculate_detunings
import ..calculate_energytrajectory
import ..calculate_normtrajectory
import ..calculate_pulses
import ..calculate_gradientsignals

function plot_trajectory(; saveas="plot_trajectory")
    isnothing(_!.work) && error("Run script before calling `plot` methods.")

    REF, FCI, FES = _!.work.system.REF, _!.work.system.FCI, _!.work.system.FES
    t = calculate_timegrid()
    E = calculate_energytrajectory()
    F = calculate_normtrajectory()

    yMAX = max(maximum(E), REF, FES)
    plot = Plots.plot(;
        xlabel= "Time (ns)",
        ylabel= "Energy (Ha)",
        xlims = axislimits(t[1], t[end]),
        ylims = axislimits(FCI, yMAX),
        legend= :topright,
        size=(1000,600),    # TODO: Fix labels out-of-bounds for non-square plot sizes.
    )

    # PLOT ENERGY, BOTH PROJECTED AND NORMALIZED
    Plots.plot!(plot, t, E, lw=3, label="Energy")
    Plots.plot!(plot, t, E./F, lw=3, label="Normalized")

    # PLOT LEAKAGE ON A SEPARATE AXIS
    twin = Plots.twinx(plot)
    Plots.plot!(twin, legend=false, ylabel="Leakage", ylims=axislimits(0,1))

    plotargs = Dict(:lw=>3, :color=>3)
    Plots.plot!(plot, [0], [2yMAX]; label="Leakage", plotargs...)   # HACK: Fix legend
    Plots.plot!(twin, t, 1 .- F, color=3, lw=3)

    # TODO: How do we get the line for E_FCI to appear above the leakage?
    #= TODO: The leakage curve is not extending all the way to t=0...
        HINT: I had to plot the hack ABOVE instead of to the LEFT of the plot,
            and that made MOST of the problem vanish. Not sure what it means.
    =#

    # PLOT DASHED LINES AT SPECIFIC ENERGIES
    Plots.hline!(plot, [REF, FCI, FES], color=:black, ls=:dash, label=false)

    !isnothing(saveas)  && Plots.savefig(plot, "$(_!.outdir)/$saveas.pdf")
    return plot
end





function plot_pulses(; f_factor=1e-2, saveas="plot_pulses")
    isnothing(_!.work) && error("Run script before calling `plot` methods.")

    t = calculate_timegrid()
    _, α, β = calculate_pulses()
    _, ϕα, ϕβ = calculate_gradientsignals()

    Ωt = plot_Ωt(t, α, β; ΩMAX=_!.setup.ΩMAX)
    ϕt = plot_ϕt(t, ϕα, ϕβ)

    # REMOVE REDUNDANT LEGENDS AND y-LABELS
    Plots.plot!(Ωt;)
    Plots.plot!(ϕt; legend=false)

    # TODO: Fix labels out-of-bounds for non-square plot sizes.
    plot = Plots.plot(Ωt, ϕt; layout = (2,1))#, size=(500,500))
    !isnothing(saveas)  && Plots.savefig(plot, "$(_!.outdir)/$saveas.pdf")
    return plot
end



function plot_trace(;
    f_calls=true, g_calls=false,
    saveas="plot_trace",
)
    isnothing(_!.work) && error("Run script before calling `plot` methods.")

    iters = _!.trace.iterations
    error = _!.trace.fn .- _!.work.system.FCI
    gnorm = _!.trace.gd

    yMAX = 2e0
    plot = Plots.plot(;
        xlabel = "Iterations",
        yscale = :log,
        ylims  = [1e-16, yMAX],
        yticks = 10.0 .^ (-16:2:0),
        legend = :bottomright,
        size=(1000,600),    # TODO: Fix labels out-of-bounds for non-square plot sizes.
    )

    Plots.plot!(plot, iters, error; lw=3, ls=:solid, label="Energy Error")
    Plots.plot!(plot, iters, gnorm; lw=3, ls=:dash, label="Gradient")

    # PLOT DASHED LINES AT CONVERGENCE CRITERIA
    convergences = [_!.setup.f_tol, _!.setup.g_tol]
    Plots.hline!(plot, convergences, color=:black, ls=:dash, label=false)

    # PLOT VERTICAL LINES AT EACH ADAPTATION
    Plots.vline!(plot, _!.trace.adaptations, color=:black, ls=:dot, label=false)

    # PLOT COUNTS ON A SEPARATE AXIS
    if f_calls || g_calls
        twin = Plots.twinx(plot)
        yMAX_ = f_calls ? last(_!.trace.f_calls) : 0
        g_calls && (yMAX_ = max(yMAX_, last(_!.trace.g_calls)))
        Plots.plot!(twin, legend=false, ylabel="Call Count", ylims=axislimits(0,yMAX_))
    end

    if f_calls
        plotargs = Dict(:lw=>2, :color=>:darkgrey, :ls=>:solid)
        Plots.plot!(plot, [0], [2yMAX]; label="f calls", plotargs...) # HACK: Fix legend
        Plots.plot!(twin, iters, _!.trace.f_calls; plotargs...)
    end

    if g_calls
        plotargs = Dict(:lw=>2, :color=>:darkgrey, :ls=>:dash)
        Plots.plot!(plot, [0], [2yMAX]; label="g calls", plotargs...) # HACK: Fix legend
        Plots.plot!(twin, iters, _!.trace.g_calls; plotargs...)
    end

    !isnothing(saveas)  && Plots.savefig(plot, "$(_!.outdir)/$saveas.pdf")
    return plot
end











function plot_Ωt(t, α, β; ΩMAX=max(maximum(abs.(α)), maximum(abs.(β))))
    yMAX = ΩMAX / 2π
    plot = Plots.plot(;
        xlabel= "Time (ns)",
        ylabel= "|Amplitude| (GHz)",
        ylims = axislimits(-yMAX, yMAX),
        legend= :topright,
    )

    # HACK: Fix legend.
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:solid, color=:black, label="α")
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:dot, color=:black, label="β")

    # PLOT AMPLITUDES
    for i in axes(α,2)
        Plots.plot!(plot, t, α[:,i]./2π, lw=3, ls=:solid, color=i, label="Drive $i")
        Plots.plot!(plot, t, β[:,i]./2π, lw=3, ls=:dot, color=i, label=false)
    end

    return plot
end

function plot_ϕt(t, ϕα, ϕβ)
    yMAX = max(maximum(abs.(ϕα)), maximum(abs.(ϕβ))) / 2π
    plot = Plots.plot(;
        xlabel= "Time (ns)",
        ylabel= "|Gradient Signal| (Ha/GHz)",
        ylims = axislimits(-yMAX, yMAX),
        legend= :topright,
    )

    # HACK: Fix legend.
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:solid, color=:black, label="ϕα")
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:dot, color=:black, label="ϕβ")

    # PLOT AMPLITUDES
    for i in axes(ϕα,2)
        Plots.plot!(plot, t, ϕα[:,i]./2π, lw=3, ls=:solid, color=i, label="Drive $i")
        Plots.plot!(plot, t, ϕβ[:,i]./2π, lw=3, ls=:dot, color=i, label=false)
    end

    return plot
end

function axislimits(MIN, MAX)
    pad = (MAX - MIN) * 0.1
    return [MIN-pad, MAX+pad]
end