import Plots

# NOTE: Imports are redundant, but provided for the sake of explicit declaration.
import ._!
import .calculate_detunings
import .calculate_trajectory
import .calculate_pulses
import .calculate_gradientsignals
import .calculate_timegrid
import .calculate_frequencygrid
import .normalized_rffts

function plot_trajectory(; saveas="plot_trajectory")
    isnothing(_!.work) && error("Run script before calling `plot` methods.")

    REF, FCI, FES = _!.work.REF, _!.work.FCI, _!.work.FES
    t = calculate_timegrid()
    E, F = calculate_trajectory()

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
    f = calculate_frequencygrid()
    Δ = calculate_detunings()

    α, β = calculate_pulses()
    α̂ = normalized_rffts(α)
    β̂ = normalized_rffts(β)

    ϕα, ϕβ = calculate_gradientsignals()
    ϕ̂α = normalized_rffts(ϕα)
    ϕ̂β = normalized_rffts(ϕβ)

    kMAX = max(
        cutoff_frequencyindex(ϕ̂α, f_factor),
        cutoff_frequencyindex(ϕ̂β, f_factor),
    )

    Ωt = plot_Ωt(t, α, β; ΩMAX=_!.setup.ΩMAX)
    Ωf = plot_Ωf(f, α̂, β̂; kMAX=kMAX, Δ=Δ)
    ϕt = plot_ϕt(t, ϕα, ϕβ)
    ϕf = plot_ϕf(f, ϕ̂α, ϕ̂β; kMAX=kMAX, Δ=Δ)

    # REMOVE REDUNDANT LEGENDS AND y-LABELS
    Plots.plot!(Ωt; legend=false)
    Plots.plot!(Ωf; ylabel="")
    Plots.plot!(ϕt; legend=false)
    Plots.plot!(ϕf; legend=false, ylabel="")

    # TODO: Fix labels out-of-bounds for non-square plot sizes.
    plot = Plots.plot(Ωt, Ωf, ϕt, ϕf; layout = (2,2), size=(1000,600))
    !isnothing(saveas)  && Plots.savefig(plot, "$(_!.outdir)/$saveas.pdf")
    return plot
end



function plot_trace(;
    gd_energy=false, gd_penalty=true, f_calls=true, g_calls=false,
    saveas="plot_trace",
)
    isnothing(_!.work) && error("Run script before calling `plot` methods.")

    iters = _!.trace.iterations
    error = _!.trace.fn_energy .- _!.work.FCI
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

    if gd_energy
        gnorm = _!.trace.gd_energy
        Plots.plot!(plot, iters, gnorm; lw=3, ls=:dot, label="Energy Gradient")
    end

    if gd_penalty
        gnorm = _!.trace.gd_penalty
        Plots.plot!(plot, iters, gnorm; lw=3, ls=:dashdot, label="Penalty Gradient")
    end

    # PLOT DASHED LINES AT CONVERGENCE CRITERIA
    Plots.hline!(plot, [_!.opt.f_tol, _!.opt.g_tol], color=:black, ls=:dash, label=false)


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

function plot_Ωf(f, α̂, β̂; kMAX=length(f), Δ=nothing)
    yMAX = max(maximum(abs.(α̂[2:end])), maximum(abs.(β̂[2:end])))
    plot = Plots.plot(;
        xlabel= "Frequency (GHz)",
        ylabel= "|Amplitude| (GHz)",
        xlims = axislimits(0, f[kMAX]),
        ylims = axislimits(0, yMAX),
        legend= :topright,
    )

    # HACK: Fix legend.
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:solid, color=:black, label="α")
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:dot,   color=:black, label="β")

    # PLOT AMPLITUDES
    k̄ = 2:kMAX  # OMIT ZERO-FREQUENCY MODE
    for i in axes(α̂,2)
        Plots.plot!(plot, f[k̄], abs.(α̂[k̄,i]), lw=3, ls=:solid, color=i, label="Drive $i")
        Plots.plot!(plot, f[k̄], abs.(β̂[k̄,i]), lw=3, ls=:dot,   color=i, label=false)
        if !isnothing(Δ)
            Plots.vline!(plot, [abs(Δ[i]/2π)], linestyle=:dash, color=i, label=false)
        end
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

function plot_ϕf(f, ϕ̂α, ϕ̂β; kMAX=length(f), Δ=nothing)
    yMAX = max(maximum(abs.(ϕ̂α[2:end])), maximum(abs.(ϕ̂β[2:end])))
    plot = Plots.plot(;
        xlabel= "Frequency (GHz)",
        ylabel= "|Gradient Signal| (Ha/GHz)",
        xlims = axislimits(0, f[kMAX]),
        ylims = axislimits(0, yMAX),
        legend= :topright,
    )

    # HACK: Fix legend.
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:solid, color=:black, label="ϕ̂α")
    Plots.plot!(plot, [0], [2yMAX], lw=3, ls=:dot, color=:black, label="ϕ̂β")

    # PLOT AMPLITUDES
    k̄ = 2:kMAX  # OMIT ZERO-FREQUENCY MODE
    for i in axes(ϕ̂α,2)
        Plots.plot!(plot, f[k̄], abs.(ϕ̂α[k̄,i]), lw=3,ls=:solid,color=i,label="Drive $i")
        Plots.plot!(plot, f[k̄], abs.(ϕ̂β[k̄,i]), lw=3, ls=:dot, color=i, label=false)
        if !isnothing(Δ)
            Plots.vline!(plot, [abs(Δ[i]/2π)], linestyle=:dash, color=i, label=false)
        end
    end

    return plot
end

function axislimits(MIN, MAX)
    pad = (MAX - MIN) * 0.1
    return [MIN-pad, MAX+pad]
end