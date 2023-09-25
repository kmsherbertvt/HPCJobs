#=

This version of AdaptiveWindows focuses on intelligent adaptations,
    passing on information about the approximate Hessian to the optimizer
    at each adaptation.

It also incorporates the revelation that my hydrogen chain Hamiltonians
    are perfectly optimizable from the *right* reference state,
    but identifying the right reference state requires the particle and spin operators.

This basically removes any good reason to randomly initialize parameters ever,
    so I am removing randomness parameters from the interface.

I am also removing the Fourier transforms from the standard plotting workflow.
I think we've pretty much learned what we needed to from them,
    and they don't do us much good anymore.

=#

module StrictAdaptiveWindows
    import Serialization: serialize, deserialize
    import Parameters: @with_kw
    import Printf: @sprintf

    import CtrlVQE
    import CtrlVQE: EvolutionType, CostFunctionType, DeviceType
    import CtrlVQE.TransmonDevices: AbstractTransmonDevice
    import CtrlVQE.TransmonDevices: resonancefrequency, drivefrequency, drivesignal

    export _!, unload, save, load, inspect, report, plot
    export iterid, adaptid, archive, unarchive

    const Float = Float64
    const as_dict(obj) = Dict(
        fld => getfield(obj, fld) for fld in fieldnames(typeof(obj))
    )

    const MATRIXPATH = "$(ENV["HOME"])/HPCJobs/matrix"
    const SYSTEMPATH = "$(ENV["HOME"])/HPCJobs/StrictAdaptiveWindows/sys"

    module MolecularSystems; include("systems.jl"); end
    import .MolecularSystems: MolecularSystem

    """ SetupVariables uniquely define a trajectory. """
    @with_kw mutable struct SetupVariables
        # SYSTEM PARAMETERS
        systemcode::String = "H2_sto-3g_singlet_1.5_P-m"
        T::Float = 10.0                     # DURATION OF PULSE

        # SIMULATION PARAMETERS
        r::Int = 20T                        # NUMBER OF TROTTER STEPS
        m::Int = 2                          # NUMBER OF LEVELS PER TRANSMON

        # HARDWARE BOUNDS
        ΩMAX::Float = 2π * 0.02 # GHz       # LARGEST PULSE AMPLITUDE
        ΔMAX::Float = 2π * 1.00 # GHz       # LARGEST PULSE DETUNING
        ΔsMIN::Float = 3.0 # ns             # SMALLEST WINDOW INTERVAL, AKIN TO BANDWIDTH

        # OPTIMIZATION TERMINATION
        f_tol::Float = 0.0
        g_tol::Float = 1e-6
        maxiter::Int = 10000
        fnRATIO::Float = 3.0                # NEVER EXCEED RATIO FN CALLS / ITERATIONS

        # ADAPTATION TERMINATION
        ΔgMIN::Float = 2g_tol
    end

    """ OptimizationVariables determine termination and update settings. """
    @with_kw mutable struct OptimizationVariables
        maxadapt::Int = 200
        update::Int = 100
    end

    """ TraceVariables store plottable data over the course of a trajectory. """
    @with_kw mutable struct TraceVariables
        iterations::Vector{Int} = Int[]
        f_calls::Vector{Int} = Int[]
        g_calls::Vector{Int} = Int[]

        adaptations::Vector{Int} = Int[]    # Stores iterations where an adaptation occurs.

        fn::Vector{Float} = Float[]
        gd::Vector{Float} = Float[]
    end

    """ StateVariables contain all the dynamic data needed to resume a trajectory. """
    @with_kw mutable struct StateVariables
        x::Vector{Float}
        Ω::Vector{Int}
        φ::Vector{Int}
        ν::Vector{Int}
        s::Vector{Vector{Float}}
    end

    """ WorkVariables contain all the fancy objects that have to be constructed. """
    @with_kw mutable struct WorkVariables
        system::MolecularSystem
        evolution::EvolutionType
        device::AbstractTransmonDevice
        energyfn::CostFunctionType
        normfn::CostFunctionType
    end

    """ Variables bundles all of the above. """
    @with_kw mutable struct Variables
        outdir::String = "job"

        setup::SetupVariables = SetupVariables()
        opt::OptimizationVariables = OptimizationVariables()
        trace::TraceVariables = TraceVariables()
        state::Union{Nothing,StateVariables} = nothing
        work::Union{Nothing,WorkVariables} = nothing

        run::Bool = true    # FLAG THAT SCRIPT SHOULD RUN LENGTHY CALCULATIONS
    end

    function unload()
        global _! = Variables()
    end
    unload()                # NOTE: This line initializes the global variable `_!`.

    function save()
        !isdir(_!.outdir) && mkdir(_!.outdir)
        serialize("$(_!.outdir)/setup", as_dict(_!.setup))
        serialize("$(_!.outdir)/opt",   as_dict(_!.opt))
        serialize("$(_!.outdir)/trace", as_dict(_!.trace))
        serialize("$(_!.outdir)/state", as_dict(_!.state))
    end

    function load(outdir; newdir=outdir, opt=true, trace=true, state=true, run=true)
        kwargs = Dict(:outdir=>newdir, :work=>nothing, :run=>run)
        kwargs[:setup] = deserialize_dict("$outdir/setup", SetupVariables)
        opt   && (kwargs[:opt]   = deserialize_dict("$outdir/opt", OptimizationVariables))
        trace && (kwargs[:trace] = deserialize_dict("$outdir/trace", TraceVariables))
        state && (kwargs[:state] = deserialize_dict("$outdir/state", StateVariables))
        global _! = Variables(; kwargs...)
        return _!
    end

    function inspect()
        notepad = """

            Setup Variables
            ---------------
        """
        for field in fieldnames(SetupVariables)
            value = getfield(_!.setup, field)
            notepad *= @sprintf "%8s: %s\n" field value
        end

        notepad *= """

            Optimization Variables
            ----------------------
        """
        for field in fieldnames(OptimizationVariables)
            value = getfield(_!.opt, field)
            notepad *= @sprintf "%8s: %s\n" field value
        end

        notepad *= """

            Ansatz
            ----------------
            Adaptations: $(length(_!.trace.adaptations))
             Parameters: $(length(_!.state.x))
        """

        notepad *= """

            Trace
            ----------------
        """
        nrecords = length(_!.trace.iterations)
        notepad *= @sprintf(
            "\n%6s %6s %6s %13s %13s\n",
            "Iter", "# f", "# ∇f", "f(x)", "∇f(x)",
        )

        nshow = min(5, nrecords)            # WE'LL SHOW ONLY THE LAST 5 RECORDS
        if nrecords > nshow                 # DEPICT THE REST WITH ELLIPSES. ^_^
            notepad *= @sprintf("%6s %6s %6s %13s %13s\n", "⋮"^5...)
        end
        for i in 1+nrecords-nshow:nrecords
            notepad *= @sprintf(
                "%6i %6i %6i %13G %13G\n",
                _!.trace.iterations[i],
                _!.trace.f_calls[i],
                _!.trace.g_calls[i],

                _!.trace.fn[i],
                _!.trace.gd[i],
            )
        end

        sysfile = "$SYSTEMPATH/$(_!.setup.systemcode)"
        if isfile(sysfile)
            system = MolecularSystems.load_system(_!.setup.systemcode)
            COR = system.REF - system.FCI
            GAP = system.FES - system.FCI
            Ex = last(_!.trace.fn)                      # LOSS FUNCTION
            εE = Ex - system.FCI                        # ENERGY ERROR
            cE = 1 - εE/COR                             # CORRELATION ENERGY RECOVERED
            gE = 1 - εE/GAP                             # "GAP" ENERGY RECOVERED
            notepad *= """

                Accuracy
                ----------------
                  Energy (Ha): $Ex
                 Energy Error: $εE
                % Corr Energy: $(cE*100)
                %  Gap Energy: $(gE*100)
            """
        end

        println(notepad)
    end

    function report()
        isnothing(_!.work) && error("`report()` needs complete setup. Try `inspect()?`")

        Δ̄ = calculate_detunings()

        notepad = "\n"
        if length(_!.trace.iterations) > 0
            system = _!.work.system
            COR = system.REF - system.FCI
            GAP = system.FES - system.FCI
            Ex = last(_!.trace.fn)                      # LOSS FUNCTION
            εE = Ex - system.FCI                        # ENERGY ERROR
            cE = 1 - εE/COR                             # CORRELATION ENERGY RECOVERED
            gE = 1 - εE/GAP                             # "GAP" ENERGY RECOVERED
            g0 = last(_!.trace.gd)                              # GRADIENT NORM

            notepad *= """

                Optimization Results
                --------------------
                  Energy (Ha): $Ex
                 Energy Error: $εE
                % Corr Energy: $(cE*100)
                %  Gap Energy: $(gE*100)
                Gradient Norm: $g0
            """
        end

        Ωtot = length(_!.state.Ω)
        Ωsat = sum((
            count(_!.state.x[_!.state.Ω] .≥    _!.setup.ΩMAX),
            count(_!.state.x[_!.state.Ω] .≤ .- _!.setup.ΩMAX),
        ))

        νtot = length(Δ̄)
        νsat = sum((
            count(Δ̄ .≥    _!.setup.ΔMAX),
            count(Δ̄ .≤ .- _!.setup.ΔMAX),
        ))

        notepad *= """

            Saturated Bounds
            ----------------
            Amplitude: $Ωsat / $Ωtot
            Frequency: $νsat / $νtot
        """

        notepad *= """

            Ansatz
            ----------------
            Adaptations: $(length(_!.trace.adaptations))
             Parameters: $(length(_!.state.x))
        """

        println(notepad)
    end

    function iterid(iter)
        return "$(lpad(iter, 8, "0"))"
    end

    function adaptid(adapt)
        return "$(lpad(adapt, 6, "0"))"
    end

    function archive(name)
        !isdir(_!.outdir) && mkdir(_!.outdir)
        serialize("$(_!.outdir)/$name", as_dict(_!.state))
    end

    function unarchive(name)
        _!.state = deserialize_dict("$(_!.outdir)/$name", StateVariables)
    end

    function deserialize_dict(file, type)
        kwargs = deserialize(file)
        isempty(kwargs) && return nothing
        return type(; kwargs...)
    end




    function calculate_timegrid()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        _, _, t̄ = CtrlVQE.trapezoidaltimegrid(_!.setup.T, _!.setup.r)
        return t̄
    end

    function calculate_detunings()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        CtrlVQE.Parameters.bind(_!.work.device, _!.state.x)

        nD = CtrlVQE.ndrives(_!.work.device)
        ω̄ = resonancefrequency.(Ref(_!.work.device), 1:nD)
        ν̄ = drivefrequency.(Ref(_!.work.device), 1:nD)
            # NOTE: USE `Ref` TO ENABLE BROADCASTING

        return ω̄ .- ν̄
    end

    function calculate_energytrajectory()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        E = Array{Float}(undef, _!.setup.r+1)   # Stores energy trajectory.

        callback = CtrlVQE.trajectory_callback(_!.work.energyfn, E)
        CtrlVQE.cost_function(_!.work.energyfn; callback=callback)(_!.state.x)

        return E
    end

    function calculate_normtrajectory()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        F = Array{Float}(undef, _!.setup.r+1)   # Stores norm trajectory (for leakage).

        callback = CtrlVQE.trajectory_callback(_!.work.normfn, F)
        CtrlVQE.cost_function(_!.work.normfn; callback=callback)(_!.state.x)

        return F
    end

    function calculate_pulses()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        CtrlVQE.Parameters.bind(_!.work.device, _!.state.x)
        t = calculate_timegrid()

        nD = CtrlVQE.ndrives(_!.work.device)
        Ω = Array{Complex{Float}}(undef, _!.setup.r+1, nD)
        for i in 1:nD
            Ω[:,i] .= drivesignal(_!.work.device, i)(t)
        end

        return Ω, real.(Ω), imag.(Ω)
    end

    function calculate_gradientsignals()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        nG = CtrlVQE.ngrades(_!.work.device)
        ϕ = Array{Float}(undef, _!.setup.r+1, nG)   # Stores gradient signals.
        CtrlVQE.grad_function(_!.work.energyfn; ϕ=ϕ)(_!.state.x)
            # NOTE: Not compatible with NormalizedEnergy.

        nD = CtrlVQE.ndrives(_!.work.device)
        ϕα = Array{Float}(undef, _!.setup.r+1, nD)
        ϕβ = Array{Float}(undef, _!.setup.r+1, nD)
        for i in 1:nD
            j = 2(i-1) + 1
            ϕα[:,i] .= ϕ[:,j  ]
            ϕβ[:,i] .= ϕ[:,j+1]
        end

        return ϕ, ϕα, ϕβ
    end





    """ Linearly interpolate roots of a series ϕ evaluated at points in t̄. """
    function identify_nodes(t̄, ϕ)
        s̄ = Float[]
        isempty(ϕ) && return s̄

        phase = sign(first(ϕ))          # TRACK IF WE ARE ABOVE OR BELOW ϕ=0 LINE
        for (i, ϕi) in enumerate(ϕ)
            sign(ϕi) == phase && continue   # SKIPS FIRST ITERATION, SO ϕ̄[i-1] IS SAFE

            # LINEAR INERPOLATION TO IDENTIFY TIME OF PHASE FLIP
            t = t̄[i]; t0 = t̄[i-1]; ϕ0 = ϕ[i-1]
            m = abs(ϕi-ϕ0)/(t-t0)           #     SLOPE OF LINE
            s = t0 + m*ϕ0                   # INTERCEPT OF LINE
            push!(s̄, s)

            phase = sign(ϕi)                # UPDATE WHICH SIDE OF ϕ=0 WE ARE ON
        end

        return s̄
    end

    """ Identify candidate window times by splitting existing windows into n pieces. """
    function nsect_windows(n, s̄0)
        Δs̄ = diff([s̄0..., _!.setup.T]) ./ n         # CANDIDATE SPLIT POINTS
        s̄ = Matrix{Float}(undef, length(s̄0), n-1)
        for i in 1:n-1
            s̄[:,i] .= s̄0 .+ (i .* Δs̄)
        end
        return reshape(s̄, :)
    end

    """ Remove candidate window times if the resulting window is too short. """
    function filter_splits!(s̄, s̄0)
        cursor = 1      # TRACKS LOCATION IN s̄0
        i = 1
        while i ≤ length(s̄)
            s = s̄[i]

            # ADVANCE CURSOR
            while cursor ≤ length(s̄0) && s̄0[cursor] < s; cursor += 1; end
            sL = get(s̄0, cursor-1, 0.0)         #  LARGEST SWITCHING TIME LESS THAN s
            sR = get(s̄0, cursor, _!.setup.T)    # SMALLEST SWITCHING TIME MORE THAN s

            # REMOVE A SPLIT IF IT CREATES A WINDOW WHICH IS TOO SMALL
            if s - sL < _!.setup.ΔsMIN || sR - s < _!.setup.ΔsMIN
                deleteat!(s̄, i)
            else
                i += 1
            end
        end
    end

    """ Duplicate a windowed signal exactly but include an extra split at s. """
    function adapt_windows(pulse, s)
        # IDENTIFY THE WINDOW WHERE s CURRENTLY FALLS
        k = findlast(s0 -> s0 < s, pulse.starttimes)
        isnothing(k) && error("`s` does not fall in any window.")
        window = deepcopy(pulse.windows[k]) # DUPLICATE PARAMETERS SO Ω(t) IS UNCHANGED

        # CREATE THE NEW PULSE
        windows = deepcopy(pulse.windows);          insert!(windows, k+1, window)
        starttimes = deepcopy(pulse.starttimes);    insert!(starttimes, k+1, s)
        return CtrlVQE.WindowedSignal(windows, starttimes)
    end

    """ Duplicate a uniform pulse as closely as possible with one more window. """
    function increment_windows(pulse, T)
        # CREATE THE NEW WINDOWS
        windows = deepcopy(pulse.windows)
        push!(windows, deepcopy(pulse.windows[end]))

        # AVERAGE PARAMETERS BETWEEN EACH WINDOW
        W = length(pulse.windows)       # ORIGINAL WINDOW COUNT
        for i in 2:W
            xL = CtrlVQE.Parameters.values(pulse.windows[i-1])
            xR = CtrlVQE.Parameters.values(pulse.windows[i])

            # WEIGHTED AVERAGE BETWEEN CONTRIBUTING WINDOWS
            x = ( (xL .* (i-1)) .+ (xR .* (W-(i-1))) ) ./ W
            CtrlVQE.Parameters.bind(windows[i], x)
        end

        # CREATE WINDOWED PULSE WITH UNIFORMLY-SPACED START TIMES
        starttimes = range(0.0, T, length(windows)+1)[1:end-1]
        return CtrlVQE.WindowedSignal(windows, starttimes)
    end






    module Plotting; include("plotting.jl"); end

    function plot(; label="plot", trajectory=true, pulses=true, trace=true)
        isnothing(_!.work) && error("Run script before calling `plot` methods.")

        trajectory && Plotting.plot_trajectory(; saveas="$(label)_trajectory")
        pulses && Plotting.plot_pulses(; saveas="$(label)_pulses")
        trace && Plotting.plot_trace(; saveas="$(label)_trace")
    end

end # module StrictAdaptiveWindows


