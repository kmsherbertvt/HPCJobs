#= Run the windowed script for increasing number of windows. =#
shift!(array, default=nothing) = isempty(array) ? default : popfirst!(array)

using UniformWindows

matrix = shift!(ARGS, "lih30")
T = parse(Float64, shift!(ARGS, "30.0"))
WMIN = parse(Int, shift!(ARGS, "1"))
WMAX = parse(Int, shift!(ARGS, "6"))

for W in WMIN:WMAX
    _!.setup = UniformWindows.SetupVariables(T=T, W=W)
    _!.outdir = "$(matrix)_$(T)_$(W)"

    include("optimize.jl")
    unload()
end
