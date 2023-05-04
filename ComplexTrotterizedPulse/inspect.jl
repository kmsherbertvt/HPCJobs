
shift!(array, default=nothing) = length(array) > 0 ? popfirst!(array) : default

jobdir = shift!(ARGS)
!isdir(jobdir) && error("$jobdir/ must exist and be pre-populated with data.")

import Serialization: deserialize
import Printf: @printf

# DISPLAY SETTINGS
if isfile("$jobdir/settings")
    settings = deserialize("$jobdir/settings")
    println("Settings:")
    for (field, value) in settings
        @printf "%20s: %s\n" field value
        # println("\t$field: $value")
    end
    println()
end

if isfile("$jobdir/x")
    x = deserialize("$jobdir/x")
    # INFER n FROM PARAMETERS (without reading in a large matrix)
    L = length(x)
    n = L % (settings[:r]+1)
    Ω = 1:L-n
    ν = 1+L-n:L

    # DISPLAY PARAMETERS
    println("Frequencies (GHz):")
    display(x[ν] ./ 2π)
    println()

    ΩRMS = √(sum(x[Ω].^2) ./ length(x[Ω]))
    Ωsat = sum((
        count(x .> settings[:ΩMAX] .- settings[:σΩ]),
        count(x .< settings[:σΩ] .- settings[:ΩMAX]),
    ))

    println("Amplitude Statistics:")
    @printf "\t            RMS Amplitude: %G GHz\n" (ΩRMS/2π)
    println("\t# of Saturated Amplitudes: $Ωsat")
    println()
end

# DISPLAY TRACE
if isfile("$jobdir/iterations")
    iterations = deserialize("$jobdir/iterations")
    trace_f_calls = deserialize("$jobdir/trace_f_calls")
    trace_g_calls = deserialize("$jobdir/trace_g_calls")

    trace_fn = deserialize("$jobdir/trace_fn")
    trace_fn_energy = deserialize("$jobdir/trace_fn_energy")
    trace_gd = deserialize("$jobdir/trace_gd")
    trace_gd_energy = deserialize("$jobdir/trace_gd_energy")
    trace_gd_penalty = deserialize("$jobdir/trace_gd_penalty")

    nrecords = length(iterations)
    nshow = min(10, nrecords)           # WE'LL SHOW ONLY THE LAST 10 RECORDS
    @printf(
        "%6s %6s %6s %13s %13s %13s %13s %13s\n",
        "Iter", "# f", "# ∇f", "f(x)", "E(x)", "∇f(x)", "∇E(x)", "∇λ(x)",
    )
    @printf("%6s %6s %6s %13s %13s %13s %13s %13s\n", "⋮"^8...)
    for i in 1+nrecords-nshow:nrecords
        @printf(
            "%6i %6i %6i %13G %13G %13G %13G %13G\n",
            iterations[i],
            trace_f_calls[i],
            trace_g_calls[i],

            trace_fn[i],
            trace_fn_energy[i],
            trace_gd[i],
            trace_gd_energy[i],
            trace_gd_penalty[i],
        )
    end
    println()
end