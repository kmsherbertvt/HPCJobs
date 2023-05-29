#= Gather data to plot, then plot it. =#

adaptjobs = ["adapteachnode", "adapteachnsect", "adaptnodes", "adaptnsect"]
W̄ = 1:6

matrix = "lih30"
T = 20.0

#= GET EXACT ENERGY =#
import NPZ: npzread
import LinearAlgebra: eigen, Hermitian
H = npzread("$(ENV["HOME"])/HPCJobs/matrix/$matrix.npy")
Λ, U = eigen(Hermitian(H))
FCI = Λ[1]

#= GATHER DATA =#
import Serialization: deserialize
# 1: parameters   2: energy error   3: iterations   4: cumulative iterations

DATA = Dict()

for job in adaptjobs
    path = "../AdaptiveWindows/jobs/$(job)_$(matrix)_$T"
    trace = deserialize("$path/trace")

    data = zeros(1+length(trace[:adaptations]), 4)
    for i in axes(data,1)
        adaptid = "$(lpad(i-1, 6, "0"))"
        iter = (i == size(data,1)) ? last(trace[:iterations]) : trace[:adaptations][i]

        state = deserialize("$path/optimized_$adaptid")

        data[i,1] = length(state[:x])
        data[i,2] = trace[:fn_energy][findfirst(trace[:iterations].==iter)] - FCI
        data[i,3] = iter - get(trace[:adaptations], i-1, 0)
        data[i,4] = iter
    end

    DATA[job] = data
end

# Uniform
prefix = "../UniformWindows/jobs/$(matrix)_$(T)_"
DATA["uniform"] = zeros(length(W̄), 4)
for i in axes(DATA["uniform"],1)
    path = "$prefix$(W̄[i])"
    trace = deserialize("$path/trace")
    state = deserialize("$path/final")

    DATA["uniform"][i,1] = length(state[:x])
    DATA["uniform"][i,2] = last(trace[:fn_energy]) - FCI
    DATA["uniform"][i,3] = last(trace[:iterations])
end
DATA["uniform"][:,4] .= cumsum(DATA["uniform"][:,3])

#= PLOT DATA =#
# data = [uniform, nodes, nsect]

labels = Dict(
    "adapteachnode" => "Adaptive (Parallel Nodes)",
    "adapteachnsect" => "Adaptive (Parallel Bisection)",
    "adaptnodes" => "Adaptive (Nodes)",
    "adaptnsect" => "Adaptive (Bisection)",
    "uniform" => "Non-Adaptive",
)

colors = Dict(
    "adapteachnode" => 1,
    "adapteachnsect" => 2,
    "adaptnodes" => 3,
    "adaptnsect" => 4,
    "uniform" => 5,
)

import Plots
yMAX = 2e0
plot = Plots.plot(;
    xlabel = "Parameters",
    ylabel = "Energy Error",
    yscale = :log,
    ylims  = [1e-16, yMAX],
    yticks = 10.0 .^ (-16:2:0),
    legend = :left,
    # size=(1000,600),    # TODO: Fix labels out-of-bounds for non-square plot sizes.
)
twin = Plots.twinx(plot)
Plots.plot!(twin, legend=false, ylabel="Iterations")

kwargs_energy  = Dict(:lw=>3, :ls=>:solid, :markershape=>:diamond)
kwargs_iter    = Dict(:lw=>1, :ls=>:dash)
kwargs_cumiter = Dict(:lw=>1, :ls=>:dot)
Plots.plot!(plot, [0], [2yMAX]; label="Energy Error", color=:black, kwargs_energy...)
Plots.plot!(plot, [0], [2yMAX]; label="Iterations", color=:black, kwargs_iter...)
Plots.plot!(plot, [0], [2yMAX]; label="Total Iterations", color=:black, kwargs_cumiter...)

for job in keys(DATA)
    Plots.plot!(plot, DATA[job][:,1], DATA[job][:,2];
            label=labels[job], color=colors[job], kwargs_energy...)
    Plots.plot!(twin, DATA[job][:,1], DATA[job][:,3];
            label=labels[job], color=colors[job], kwargs_iter...)
    Plots.plot!(twin, DATA[job][:,1], DATA[job][:,4];
            label=labels[job], color=colors[job], kwargs_cumiter...)
end

Plots.savefig("$(matrix)_$(T).pdf")
Plots.gui()

Plots.plot!(plot, title="$(T)ns Pulse Duration")
Plots.savefig("$(matrix)_$(T).png")