module BFGSTrotterized
    import Serialization: serialize, deserialize
    import Parameters: @with_kw
    import Printf: @sprintf

    import CtrlVQE
    import CtrlVQE: EvolutionType, CostFunctionType, DeviceType
    import CtrlVQE.TransmonDevices: AbstractTransmonDevice
    import CtrlVQE.TransmonDevices: resonancefrequency, drivefrequency, drivesignal

    export _!, unload, save, load, inspect, report, plot
    export iterid, archive, unarchive

    const Float = Float64
    const as_dict(obj) = Dict(
        fld => getfield(obj, fld) for fld in fieldnames(typeof(obj))
    )

    const MATRIXPATH = "$(ENV["HOME"])/HPCJobs/matrix"
    const SYSTEMPATH = "$(ENV["HOME"])/HPCJobs/BFGSTrotterized/sys"

    module MolecularSystems; include("systems.jl"); end
    import .MolecularSystems: MolecularSystem

    """ SetupVariables uniquely define a trajectory. """
    @with_kw mutable struct SetupVariables
        # SYSTEM PARAMETERS
        systemcode::String = "H2_sto-3g_singlet_1.5_P-m"
        T::Float = 10.0                     # DURATION OF PULSE

        # SIMULATION PARAMETERS
        r::Int = round(Int,20T)             # NUMBER OF TROTTER STEPS
        m::Int = 2                          # NUMBER OF LEVELS PER TRANSMON

        # HARDWARE BOUNDS
        ΩMAX::Float = 2π * 0.02 # GHz       # LARGEST PULSE AMPLITUDE
        ΔMAX::Float = 2π * 1.00 # GHz       # LARGEST PULSE DETUNING

        # SOFTWARE BOUNDS
        λΩ::Float = 1.0 # Ha                # STRENGTH OF PENALTY TERMS
        λΔ::Float = 1.0 # Ha
        σΩ::Float = ΩMAX                    # SMOOTHNESS OF PENALTY TERMS
        σΔ::Float = ΔMAX
    end

    """ OptimizationVariables determine termination and update settings. """
    @with_kw mutable struct OptimizationVariables
        f_tol::Float = 0.0
        g_tol::Float = 1e-6
        maxiter::Int = 10000
        fnRATIO::Float = 10.0               # NEVER EXCEED RATIO FN CALLS / ITERATIONS
        update::Int = 500
    end

    """ TraceVariables store plottable data over the course of a trajectory. """
    @with_kw mutable struct TraceVariables
        iterations::Vector{Int} = Int[]
        f_calls::Vector{Int} = Int[]
        g_calls::Vector{Int} = Int[]

        fn::Vector{Float} = Float[]
        gd::Vector{Float} = Float[]
        energyfn::Vector{Float} = Float[]
        boundsfn::Vector{Float} = Float[]
        energygd::Vector{Float} = Float[]
        boundsgd::Vector{Float} = Float[]
    end

    """ StateVariables contain all the dynamic data needed to resume a trajectory. """
    @with_kw mutable struct StateVariables
        x::Vector{Float}
        Ω::Vector{Int}
        φ::Vector{Int}
        ν::Vector{Int}
        s::Vector{Vector{Float}}
        Hk::Matrix{Float}           # INVERSE HESSIAN APPROXIMATION
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
             Parameters: $(length(_!.state.x))
        """

        nrecords = length(_!.trace.iterations)
        notepad *= @sprintf(
            "\n%6s %6s %6s %13s %13s %13s %13s %13s %13s\n",
            "Iter", "# f", "# ∇f", "f(x)", "∇f(x)", "E(x)", "∇E(x)", "∇E(x)", "∇λ(x)",
        )

        nshow = min(5, nrecords)            # WE'LL SHOW ONLY THE LAST 5 RECORDS
        if nrecords > nshow                 # DEPICT THE REST WITH ELLIPSES. ^_^
            notepad *= @sprintf("%6s %6s %6s %13s %13s %13s %13s %13s %13s\n", "⋮"^9...)
        end

        for i in 1+nrecords-nshow:nrecords
            notepad *= @sprintf(
                "%6i %6i %6i %13G %13G %13G %13G %13G %13G\n",
                _!.trace.iterations[i],
                _!.trace.f_calls[i],
                _!.trace.g_calls[i],

                _!.trace.fn[i],
                _!.trace.gd[i],
                _!.trace.energyfn[i],
                _!.trace.energygd[i],
                _!.trace.boundsfn[i],
                _!.trace.boundsgd[i],
            )
        end

        sysfile = "$SYSTEMPATH/$(_!.setup.systemcode)"
        if isfile(sysfile)
            system = MolecularSystems.load_system(_!.setup.systemcode)
            COR = system.REF - system.FCI
            GAP = system.FES - system.FCI
            Ex = last(_!.trace.energyfn)                # LOSS FUNCTION
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
            Ex = last(_!.trace.energyfn)                # BEST ENERGY
            εE = Ex - system.FCI                        # ENERGY ERROR
            cE = 1 - εE/COR                             # CORRELATION ENERGY RECOVERED
            gE = 1 - εE/GAP                             # "GAP" ENERGY RECOVERED

            f0 = last(_!.trace.fn)                      # LOSS FUNCTION
            λx = f0 - Ex                                # PENALTY CONTRIBUTION
            g0 = last(_!.trace.gd)                      # GRADIENT NORM

            notepad *= """

                Optimization Results
                --------------------
                  Energy (Ha): $Ex
                 Energy Error: $εE
                % Corr Energy: $(cE*100)
                %  Gap Energy: $(gE*100)

                Loss Function: $f0
                 from  Energy: $Ex
                 from Penalty: $λx
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
             Parameters: $(length(_!.state.x))
        """

        println(notepad)
    end

    function iterid(iter)
        return "$(lpad(iter, 8, "0"))"
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




    module Plotting; include("plotting.jl"); end

    function plot(; label="plot", trajectory=true, pulses=true, trace=true)
        isnothing(_!.work) && error("Run script before calling `plot` methods.")

        trajectory && Plotting.plot_trajectory(; saveas="$(label)_trajectory")
        pulses && Plotting.plot_pulses(; saveas="$(label)_pulses")
        trace && Plotting.plot_trace(; saveas="$(label)_trace")
    end

end # module BFGSTrotterized
