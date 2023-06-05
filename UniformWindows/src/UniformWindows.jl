module UniformWindows
    import Serialization: serialize, deserialize
    import Parameters: @with_kw
    import Printf: @sprintf

    import FFTW: rfft, rfftfreq

    import CtrlVQE
    import CtrlVQE.CostFunctions: AbstractCostFunction, AbstractGradientFunction
    import CtrlVQE.Devices.TransmonDevices: AbstractTransmonDevice
    import CtrlVQE.Devices.TransmonDevices: resonancefrequency,drivefrequency,drivesignal

    export _!, unload, save, load, inspect, report, iterid, archive, unarchive, plot

    const Float = Float64

    @with_kw mutable struct SetupVariables
        matrix::String = "lih30"
        T::Float = 20.0 # ns
        W::Int = 10

        r::Int = 20T
        m::Int = 2

        seed::Int = 9999
        init_Ω::Float = 0.0 # 2π GHz
        init_φ::Float = 0.0
        init_Δ::Float = 0.0 # 2π GHz

        ΩMAX::Float = 2π * 0.02 # GHz
        λΩ::Float = 1.0 # Ha
        σΩ::Float = ΩMAX

        ΔMAX::Float = 2π * 1.00 # 2π GHz
        λΔ::Float = 1.0 # Ha
        σΔ::Float = ΔMAX
    end

    @with_kw mutable struct OptimizationVariables
        f_tol::Float = 0.0
        g_tol::Float = 1e-6
        maxiter::Int = 10000
        fnRATIO::Float = 3.0
        update::Int = 100
    end

    @with_kw mutable struct TraceVariables
        iterations::Vector{Int} = Int[]
        f_calls::Vector{Int} = Int[]
        g_calls::Vector{Int} = Int[]

        fn::Vector{Float} = Float[]
        fn_energy::Vector{Float} = Float[]
        gd::Vector{Float} = Float[]
        gd_energy::Vector{Float} = Float[]
        gd_penalty::Vector{Float} = Float[]
    end

    @with_kw mutable struct StateVariables
        x::Vector{Float}
        Ω::Vector{Int}
        φ::Vector{Int}
        ν::Vector{Int}
    end

    @with_kw mutable struct WorkVariables
        REF::Float
        FCI::Float
        FES::Float

        ψ_FCI::Vector{Complex{Float}}

        device::AbstractTransmonDevice
        fn_energy::AbstractCostFunction
        gd_energy::AbstractGradientFunction
        fn_norm::AbstractCostFunction
    end

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
        serialize_dict("$(_!.outdir)/setup", _!.setup)
        serialize_dict("$(_!.outdir)/opt",   _!.opt)
        serialize_dict("$(_!.outdir)/trace", _!.trace)
        serialize_dict("$(_!.outdir)/state", _!.state)
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

        nrecords = length(_!.trace.iterations)
        notepad *= @sprintf(
            "\n%6s %6s %6s %13s %13s %13s %13s %13s\n",
            "Iter", "# f", "# ∇f", "f(x)", "E(x)", "∇f(x)", "∇E(x)", "∇λ(x)",
        )

        nshow = min(5, nrecords)            # WE'LL SHOW ONLY THE LAST 5 RECORDS
        if nrecords > nshow                 # DEPICT THE REST WITH ELLIPSES. ^_^
            notepad *= @sprintf("%6s %6s %6s %13s %13s %13s %13s %13s\n", "⋮"^8...)
        end

        for i in 1+nrecords-nshow:nrecords
            notepad *= @sprintf(
                "%6i %6i %6i %13G %13G %13G %13G %13G\n",
                _!.trace.iterations[i],
                _!.trace.f_calls[i],
                _!.trace.g_calls[i],

                _!.trace.fn[i],
                _!.trace.fn_energy[i],
                _!.trace.gd[i],
                _!.trace.gd_energy[i],
                _!.trace.gd_penalty[i],
            )
        end

        println(notepad)
    end

    function report()
        isnothing(_!.work) && error("`report()` needs complete setup. Try `inspect()?`")

        Δ̄ = calculate_detunings()

        notepad = "\n"
        if length(_!.trace.iterations) > 0
            f0 = last(_!.trace.fn)              # LOSS FUNCTION
            E0 = last(_!.trace.fn_energy)       # CURRENT ENERGY
            λ0 = f0 - E0                        # PENALTY CONTRIBUTION
            εE = E0 - _!.work.FCI               # ENERGY ERROR
            cE = 1-εE/(_!.work.REF-_!.work.FCI) # CORRELATION ENERGY
            g0 = last(_!.trace.gd)              # GRADIENT NORM

            notepad *= """

                Optimization Results
                --------------------
                 Energy Error: $εE
                % Corr Energy: $(cE*100)

                Loss Function: $f0
                 from  Energy: $E0
                 from Penalty: $λ0

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

        println(notepad)
    end

    function iterid(iter)
        return "state_$(lpad(iter, 6, "0"))"
    end

    function archive(name)
        !isdir(_!.outdir) && mkdir(_!.outdir)
        serialize_dict("$(_!.outdir)/$name", _!.state)
    end

    function unarchive(name)
        _!.state = deserialize_dict("$(_!.outdir)/$name", StateVariables)
    end






    function serialize_dict(file, obj)
        kwargs = Dict(field => getfield(obj, field) for field in fieldnames(typeof(obj)))
        serialize(file, kwargs)
    end

    function deserialize_dict(file, type, null=false)
        kwargs = deserialize(file)
        isempty(kwargs) && return nothing
        return type(; kwargs...)
    end


    function calculate_detunings()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        nD = CtrlVQE.ndrives(_!.work.device)
        ω̄ = resonancefrequency.(Ref(_!.work.device), 1:nD)
        ν̄ = drivefrequency.(Ref(_!.work.device), 1:nD)

        return ω̄ .- ν̄
    end

    function calculate_trajectory()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")


        E = Array{Float}(undef, _!.setup.r+1)   # Stores energy trajectory.
        F = Array{Float}(undef, _!.setup.r+1)   # Stores norm trajectory (for leakage).

        CtrlVQE.Parameters.bind(_!.work.device, _!.state.x)
        # CtrlVQE.evolve(
        #     # All this part is just boiler-plate to make an evolution happen.
        #     _!.work.fn_energy.algorithm,
        #     _!.work.fn_energy.device,
        #     _!.work.fn_energy.basis,
        #     _!.work.fn_energy.T,
        #     _!.work.fn_energy.ψ0;
        #     result=_!.work.fn_energy.ψ,
        # TODO: ^- This doesn't work. Something's wrong with callback in `evolve`?
        CtrlVQE.gradientsignals(
            # All this part is just boiler-plate to make an evolution happen.
            _!.work.gd_energy.f.device,
            _!.work.gd_energy.f.basis,
            _!.work.gd_energy.f.T,
            _!.work.gd_energy.f.ψ0,
            _!.work.gd_energy.r,
            _!.work.gd_energy.f.OT;
            result=_!.work.gd_energy.ϕ̄,
            evolution=_!.work.gd_energy.f.algorithm,

            # The important part is here: fill out energy and norm arrays at each step.
            callback=(i, t, ψ) -> (
                E[i] = CtrlVQE.evaluate(_!.work.fn_energy, ψ, t);
                F[i] = CtrlVQE.evaluate(_!.work.fn_norm,   ψ, t);
            ),
        )

        return E, F
    end

    function calculate_pulses()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        CtrlVQE.Parameters.bind(_!.work.device, _!.state.x)
        t = calculate_timegrid()

        nD = CtrlVQE.ndrives(_!.work.device)
        α = Array{Float}(undef, _!.setup.r+1, nD)
        β = Array{Float}(undef, _!.setup.r+1, nD)
        for i in 1:nD
            Ωt = drivesignal(_!.work.device, i)(t)
            α[:,i] = real.(Ωt)
            β[:,i] = imag.(Ωt)
        end

        return α, β
    end

    function calculate_gradientsignals()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        # RUN ENERGY GRADIENT FUNCTION TO FILL OUT GRADIENT SIGNALS
        CtrlVQE.Parameters.bind(_!.work.device, _!.state.x)
        CtrlVQE.gradientsignals(
            _!.work.gd_energy.f.device,
            _!.work.gd_energy.f.basis,
            _!.work.gd_energy.f.T,
            _!.work.gd_energy.f.ψ0,
            _!.work.gd_energy.r,
            _!.work.gd_energy.f.OT;
            result=_!.work.gd_energy.ϕ̄,
            evolution=_!.work.gd_energy.f.algorithm,
        )
        ϕ = _!.work.gd_energy.ϕ̄

        nD = CtrlVQE.ndrives(_!.work.device)
        ϕα = Array{Float}(undef, _!.setup.r+1, nD)
        ϕβ = Array{Float}(undef, _!.setup.r+1, nD)
        for i in 1:nD
            j = 2(i-1) + 1
            ϕα[:,i] .= ϕ[:,j  ]
            ϕβ[:,i] .= ϕ[:,j+1]
        end

        return ϕα, ϕβ
    end

    function calculate_timegrid()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        _, _, t̄ = CtrlVQE.trapezoidaltimegrid(_!.setup.T, _!.setup.r)
        return t̄
    end

    function calculate_frequencygrid()
        isnothing(_!.work) && error("Run script before calling `calculate` methods.")

        τ, _, _ = CtrlVQE.trapezoidaltimegrid(_!.setup.T, _!.setup.r)
        return rfftfreq(_!.setup.r, 1/τ)
    end

    function normalized_rffts(fs)   # `fs` is a matrix, with each column a different f(t).
        fs = fs[2:end, :]           # Skip t=0 so time grid is uniform.
        f̂s = rfft(fs, 1)            # `rfft` omits redundant negative frequencies.
        f̂s ./= size(fs,1)           # Standard FFT normalization.
        f̂s[2:end] .*= 2             # Account for mass "lost" by omitting negative freqs.
        return f̂s
    end

    function cutoff_frequencyindex(f̂s, f_factor)
        cutoff = maximum(abs.(f̂s[2:end])) * f_factor
        index = 1
        default(x,y) = isnothing(x) ? y : x
        for i in axes(f̂s, 2)
            index = max(index, default(findlast(abs.(f̂s[:,i]) .> cutoff), 0))
        end
        return index
    end


    function plot(; label="plot", trajectory=true, pulses=true, trace=true)
        isnothing(_!.work) && error("Run script before calling `plot` methods.")

        trajectory && plot_trajectory(; saveas="$(label)_trajectory")
        pulses && plot_pulses(; saveas="$(label)_pulses")
        trace && plot_trace(; saveas="$(label)_trace")
    end

    include("plotting.jl")          # Provides methods to create individual subplots.

end # module UniformWindows


