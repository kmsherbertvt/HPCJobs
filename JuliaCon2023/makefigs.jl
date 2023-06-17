#= Script to load the data and generate the plots we want to use. =#

jobwithfreq = "optimizewithfreq_lih30_100.0_50"
jobsansfreq = "optimize_lih30_100.0_50"

ΩMAX = 2π * 0.002
tMIN = 5.0
tMAX = 15.0
kMAX = 51

using JuliaCon2023
import Plots

load(jobwithfreq; run=false)
try; include("$jobwithfreq/script.jl"); catch; end
Δi = _!.work.device.ω̄ - xi[ν]

Ωt, ϕt, ϕf = JuliaCon2023.plot_pulse_figures(iterid(0); tMIN=tMIN, tMAX=tMAX, ΩMAX=ΩMAX)
Plots.savefig(Ωt, "fig/init_Ωt.pdf")
Plots.savefig(ϕt, "fig/init_ϕt.pdf")

Ωt, ϕt, ϕf = JuliaCon2023.plot_pulse_figures(iterid(1); tMIN=tMIN, tMAX=tMAX, ΩMAX=ΩMAX)
Plots.savefig(Ωt, "fig/step_Ωt.pdf")
Plots.savefig(ϕt, "fig/step_ϕt.pdf")

lastiter = _!.trace.iterations[end]
println("With Frequencies: $lastiter Iterations")
Ωt, ϕt, ϕf = JuliaCon2023.plot_pulse_figures(iterid(lastiter);
    kMAX=51, tMIN=tMIN, tMAX=tMAX)
Plots.savefig(Ωt, "fig/final_Ωt.pdf")
Plots.savefig(ϕt, "fig/final_ϕt.pdf")
Plots.savefig(ϕf, "fig/final_ϕf_sans_Δ.pdf")
Ωt, ϕt, ϕf = JuliaCon2023.plot_pulse_figures(iterid(lastiter); kMAX=51, plot_Δ=true)
Plots.savefig(Ωt, "fig/final_Ωt_full.pdf")
Plots.savefig(ϕf, "fig/final_ϕf_with_Δ.pdf")
for i in 1:n
    Plots.vline!(ϕf, [abs(Δi[i]/2π)], linestyle=:dot, color=i, label=false)
end
Plots.savefig(ϕf, "fig/final_ϕf_with_Δ_and_Δi.pdf")


load(jobsansfreq; run=false)
try; include("$jobsansfreq/script.jl"); catch; end
lastiter = _!.trace.iterations[end]
println("Sans Frequencies: $lastiter Iterations")
Ωt, ϕt, ϕf = JuliaCon2023.plot_pulse_figures(iterid(lastiter))
Plots.savefig(Ωt, "fig/resonant_Ωt_full.pdf")