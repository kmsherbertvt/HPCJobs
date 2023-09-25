#= resonant_uniform

Resonant: Fix pulse frequencies on resonance.
n-sects Parallel: Split each window into a rational division at each adaption.
Complex: Use complex pulses with inscribed amplitude bound.
Safe: If a parameter is adapted, reset its part of the inverse Hessian to identity.
Smooth: Use smooth bounds for amplitude and drive frequency.

=#

import BFGSAdaptiveWindows as JOB
import BFGSAdaptiveWindows: _!

import CtrlVQE

import LinearAlgebra: eigen, Hermitian, norm, diagm

using Serialization: serialize, deserialize
using NPZ: npzread
import Optim, LineSearches

##########################################################################################
#= INITIALIZE SETUP VARIABLES =#

# LOAD THE MOLECULAR SYSTEM (INCLUDING H, REFERENCE STATE, EXACT GROUNDSTATE)
system = JOB.MolecularSystems.load_system(_!.setup.systemcode)

# PREP EVOLUTION ALGORITHM
evolution = CtrlVQE.Toggle(_!.setup.r)

# PREP DEVICE - THE PULSE IS A PLACEHOLDER UNTIL CURRENT STATE IS LOADED
protopulse = CtrlVQE.ComplexConstant(zero(JOB.Float), zero(JOB.Float))
device = CtrlVQE.Systematic(
    CtrlVQE.FixedFrequencyTransmonDevice,
    system.n,
    CtrlVQE.UniformWindowed(protopulse, _!.setup.T, 1);
    m = _!.setup.m,
)

xi = CtrlVQE.Parameters.values(device)  # FOR CONVENIENCE

# PREP ENERGY FUNCTIONS (note state changes are captured in device)
basis = CtrlVQE.OCCUPATION
frame = CtrlVQE.STATIC
ψ0 = CtrlVQE.QubitOperators.project(system.ψ_REF, device)
O0 = CtrlVQE.QubitOperators.project(system.H, device)

energyfn = CtrlVQE.ProjectedEnergy(evolution, device, basis, frame, _!.setup.T, ψ0, O0)
normfn = CtrlVQE.Normalization(evolution, device, basis, _!.setup.T, ψ0)

##########################################################################################
#= INITIALIZE STATE VARIABLES, OR UPDATE WORK VARIABLES WITH LOADED STATE VARIABLES =#

if isnothing(_!.state)
    L = length(xi)
    _!.state = JOB.StateVariables(
        x = copy(xi),
        Ω = collect(1:L), φ = [], ν = [],                   # COMPLEX RESONANT
        s = [device.Ω̄[i].starttimes for i in 1:CtrlVQE.ndrives(device)],
        Hk = diagm(ones(JOB.Float, L)),
    )
else
    cursor = 1
    for i in 1:CtrlVQE.ndrives(device)
        windows = CtrlVQE.ComplexConstant{JOB.Float}[]
        for (k, s) in enumerate(_!.state.s[i])
            push!(windows, CtrlVQE.ComplexConstant(
                _!.state.x[cursor],
                _!.state.x[cursor+1],
            ))
            global cursor += 2
        end
        device.Ω̄[i] = CtrlVQE.WindowedSignal(windows, _!.state.s[i])
    end
    CtrlVQE.Parameters.bind(device, _!.state.x)
end

##########################################################################################
#= ASSIGN WORK VARIABLES =#

_!.work = JOB.WorkVariables(
    system = system,
    evolution = evolution,
    device = device,
    energyfn = energyfn,
    normfn = normfn,
)


##########################################################################################
#= DEFINE OPTIMIZATION PROTOCOL =#

linesearch = LineSearches.MoreThuente()

function bounds(L, ω̄)
    Ω = collect(1:L); φ = []; ν = []                    # COMPLEX RESONANT

    λ  = zeros(L);
        λ[Ω] .= _!.setup.λΩ;
        λ[ν] .= _!.setup.λΔ
    μR = zeros(L); !isempty(ν) && (μR[ν] .= ω̄)
        μR[Ω] .+= _!.setup.ΩMAX ./ √2;
        μR[ν] .+= _!.setup.ΔMAX
    μL = zeros(L); !isempty(ν) && (μL[ν] .= ω̄)
        μL[Ω] .-= _!.setup.ΩMAX ./ √2;
        μL[ν] .-= _!.setup.ΔMAX
    σ  = zeros(L);
        σ[Ω] .= _!.setup.σΩ;
        σ[ν] .= _!.setup.σΔ
    return CtrlVQE.SmoothBound(λ, μR, μL, σ)
end

iteration = isempty(_!.trace.iterations) ? Ref(0) : Ref(last(_!.trace.iterations))
function callback(history)
    iteration[] += 1

    optstate = last(history)        # EXTRACT THE CURRENT OPTIMIZATION STATE...
    empty!(history)                 # ...THEN RELEASE THE REST OF THE HISTORY FROM MEMORY
    push!(history, optstate)        # ...BUT LET IT ALWAYS KEEP THE LAST ONE...

    # PRINT ITERATION LINE
    println("$(optstate.iteration)\t$(optstate.value)\t$(optstate.g_norm)")

    # UPDATE TRACES
    push!(_!.trace.iterations, iteration[])
    push!(_!.trace.f_calls,    lossfn.f_counter[])
    push!(_!.trace.g_calls,    lossfn.g_counter[])

    push!(_!.trace.fn,         optstate.value)
    push!(_!.trace.gd,         optstate.g_norm)

    push!(_!.trace.energyfn,   lossfn.values[1])
    push!(_!.trace.energygd,   lossfn.gnorms[1])
    push!(_!.trace.boundsfn,   lossfn.values[2])
    push!(_!.trace.boundsgd,   lossfn.gnorms[2])

    # UPDATE DATA
    _!.state.x  .= optstate.metadata["x"]
    _!.state.Hk .= optstate.metadata["~inv(H)"]
    JOB.save()
    if iteration[] % _!.opt.update == 0
        JOB.report()
        JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
    end

    # CHECK FOR SPECIAL TERMINATION CONDITIONS
    iteration[] > 10 && lossfn.f_counter[] > _!.setup.fnRATIO * iteration[] && (
        println("Linesearch excessively difficult");
        return true;
    )
    return false
end

options = Optim.Options(
    f_tol = _!.setup.f_tol,
    g_tol = _!.setup.g_tol,
    iterations = _!.setup.maxiter,
    callback = callback,
    store_trace = true,             # AFAIK, THIS IS THE ONLY WAY TO ACCESS Hk
    extended_trace = true,
)

dfdf = Ref{Any}(nothing)
function do_optimization()
    boundsfn = bounds(length(_!.state.x), _!.work.device.ω̄)
    global lossfn = CtrlVQE.CompositeCostFunction(_!.work.energyfn, boundsfn)

    # KEEP TRACK OF THE NUMBER OF FUNCTION CALLS
    lossfn.f_counter[] = isempty(_!.trace.f_calls) ? 0 : last(_!.trace.f_calls)
    lossfn.g_counter[] = isempty(_!.trace.g_calls) ? 0 : last(_!.trace.g_calls)


    f  = CtrlVQE.cost_function(lossfn)
    g! = CtrlVQE.grad_function_inplace(lossfn)

    Hk = _!.state.Hk
    optimizer = Optim.BFGS(
        linesearch=linesearch,
        initial_invH=(x->Hk),
    )

    optimization = Optim.optimize(
        f,
        g!,
        _!.state.x,
        optimizer,
        options,
    )
    println(optimization)

    _!.state.x  .= Optim.minimizer(optimization)
    _!.state.Hk .= last(Optim.trace(optimization)).metadata["~inv(H)"]
    return Optim.converged(optimization)
end

##########################################################################################
#= DEFINE ADAPTATION PROTOCOL =#

function map_indices(sBEST, ΔgBEST, sCURR)
    # sBEST is a vector of candidate split times, one for each pulse
    #   If the split time is zero, ignore that pulse.
    # sCURR is a vector of vectors with the current split time
    # We are constructing the `imap` to be used in `adapt_Hk`.

    println("--- ADAPTING Hk ---")
    println(sBEST)
    println(ΔgBEST)
    println(sCURR)

    n_pperw = CtrlVQE.Parameters.count(protopulse)      # PARAMETERS PER WINDOW
    n_pulse = length(sCURR)
    n_windw = [length(sCURR[i]) for i in 1:n_pulse]     # NUMBER OF WINDOWS FOR EACH PULSE
    n_param = n_pperw * sum(n_windw)                    # RESONANT (otherwise add n_pulse)
    imap = collect(1:n_param)

    # IDENTIFY INDICES WHICH WILL BE SPLIT
    i_split = Int[]
    n_cumwindow  = cumsum([0, n_windw...])  # ith ELEMENT GIVES # OF WINDOWS BEFORE PULSE i

    for i in 1:CtrlVQE.ndrives(device)    # TETRIS LOOP (PARALLEL)
    # for i in argmax(ΔgBEST)                 # OPTIMAL "LOOP" (only one iteration)
        sBEST[i] == 0 && continue   # NO s GAVE SUFFICIENT IMPROVEMENT FOR THIS PULSE

        k = findlast(s_ -> s_ < sBEST[i], sCURR[i]) # LOCATE INDEX OF WINDOW TO BE SPLIT
        for j in 1:n_pperw                  # ACCOUNT FOR EACH PARAMETER IN TARGET WINDOW
            push!(i_split, (n_cumwindow[i] + (k-1))*n_pperw + j)
        end
    end

    println(i_split)

    # RESET SPLIT PARAMETERS, AND INSERT NEW PARAMETERS
    for i in i_split
        k = findfirst(imap .== i)
        @assert !isnothing(k)               # imap HAS EVERY POSSIBLE i BY DEFINITION

        imap[k] = 0                         # SPLIT PARAMETER SHOULD BE RESET (SAFE)
        insert!(imap, k, 0)                 # EACH SPLIT PARAMETER GETS ANOTHER, ADJACENT

        #= NOTE: This itches.

        The inserted parameter is not really adjacent to the split parameter.
        All the inserted parameters for a *given window*
            are adjacent to the split parameters for the original window,
            but all the *other parameters for the original window* are in between,
            say, the first split parameter and its inserted parameter.
        But since *all* of these are set to 0, we can shortcut it like this.
        But this will be harder for a less SAFE choice of Hessian update.

        =#
    end

    return imap
end

function adapt_Hk(Hk, imap)
    # New Hk at index i corresponds to old Hk at index imap[i].
    # If imap[i] = 0, it is a brand new parameter.
    # For "safe" adaptation, treat any split parameter as brand new.

    new_Hk = diagm(ones(eltype(Hk), length(imap)))
    for ix in CartesianIndices(Hk); i,j = Tuple(ix)
        i0 = imap[i]; i0 == 0 && continue
        j0 = imap[j]; j0 == 0 && continue
        new_Hk[i,j] = Hk[i0,j0]
    end
    return new_Hk
end

function do_adaptation()
    CtrlVQE.Parameters.bind(device, _!.state.x)
    ϕ, _, _ = JOB.calculate_gradientsignals()
    _, τ̄, t̄ = CtrlVQE.trapezoidaltimegrid(_!.setup.T, _!.setup.r)

    # CALCULATE OPTIMAL SWITCHING TIMES FOR EACH PULSE
    ΔgBEST = zeros(CtrlVQE.ndrives(device))  # TRACK LARGEST Δg FOR EACH PULSE
    sBEST  = zeros(CtrlVQE.ndrives(device))  # TRACK s GIVING LARGEST Δg " "

    maxnsect = floor(Int, _!.setup.T / _!.setup.ΔsMIN)
    for nsect in 2:maxnsect
        i_tosplit = findall(s -> s == 0, sBEST) # WHICH DRIVES STILL NEED A SPLIT?
        isempty(i_tosplit) && break             # STOP ONCE ALL DRIVES HAVE A SPLIT

        for j in 1:CtrlVQE.ngrades(device)
            i = (j-1)÷2 + 1                     # IDENTIFY DRIVE
            !(i in i_tosplit) && continue       # SKIP DRIVE SPLITTABLE WITH SMALLER nsect

            pulse = device.Ω̄[i]
            s̄ = JOB.nsect_windows(nsect, pulse.starttimes)
            JOB.filter_splits!(s̄, pulse.starttimes)

            for s in s̄
                # CONSTRUCT TRIAL DEVICE
                candidate = deepcopy(device)
                candidate.Ω̄[i] = JOB.adapt_windows(device.Ω̄[i], s)
                g = CtrlVQE.Devices.gradient(candidate, τ̄, t̄, ϕ)

                # TAKE BOUNDS INTO ACCOUNT WHEN SELECTING OPTIMAL SPLIT
                candidatebounds = bounds(CtrlVQE.Parameters.count(candidate), candidate.ω̄)
                g .+= CtrlVQE.grad_function(candidatebounds)(
                    CtrlVQE.Parameters.values(candidate)
                )

                # CALCULATE CHANGE IN GRADIENT
                Δg = maximum(abs.(g))
                Δg < _!.setup.ΔgMIN && continue     # ABANDON TRIAL IF IT IS BELOW THRESHOLD

                # UPDATE OPTIMAL SWITCHING TIMES
                if Δg > ΔgBEST[i]
                    ΔgBEST[i] = Δg
                    sBEST[i]  = s
                end
            end
        end
    end

    # CONSTRUCT AN UPDATED HESSIAN
    new_Hk = adapt_Hk(_!.state.Hk, map_indices(sBEST, ΔgBEST, _!.state.s))  # SAFE

    # UPDATE THE DEVICE PULSES
    adapted = false
    for i in 1:CtrlVQE.ndrives(device)    # TETRIS LOOP (PARALLEL)
    # for i in argmax(ΔgBEST)                 # OPTIMAL "LOOP" (only one iteration)
        sBEST[i] == 0 && continue   # NO s GAVE SUFFICIENT IMPROVEMENT FOR THIS PULSE
        adapted = true              # *SOME* CHANGE WAS MADE IN THIS ADAPTATION

        device.Ω̄[i] = JOB.adapt_windows(device.Ω̄[i], sBEST[i])
    end
    !adapted && return false

    # UPDATE THE STATE
    x = CtrlVQE.Parameters.values(device)
    L = length(x)
    _!.state = JOB.StateVariables(
        x = x,
        Ω = collect(1:L), φ = [], ν = [],                   # COMPLEX RESONANT
        s = [device.Ω̄[i].starttimes for i in 1:CtrlVQE.ndrives(device)],
        Hk = new_Hk,
    )

    # UPDATE THE TRACE
    iter = length(_!.trace.iterations) > 0 ? last(_!.trace.iterations) : 0
    push!(_!.trace.adaptations, iter)

    return true
end



##########################################################################################
#= RUN OPTIMIZATION =#

!_!.run && error("Setup Finished")

JOB.save()
JOB.report()

name = "init"
JOB.archive(name)
JOB.plot(; label=name, trace=false)

converged = false
equilibrate_round = true
while converged || equilibrate_round
    # ARCHIVE OPTIMIZED STATE (skips first iteration)
    if converged
        JOB.save()
        JOB.report()

        global name = "optimized_$(JOB.adaptid(length(_!.trace.adaptations)))"
        JOB.archive(name)
        JOB.plot(; label=name, trace=false)
        JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
    end

    # EXCEEDED BOUNDS OF OPTIMIZATION VARIABLES (user can just change them)
    if length(_!.trace.adaptations) ≥ _!.opt.maxadapt
        println("Maximum adaptations reached.")
        break
    end

    # ADAPT ANSATZ AND CHECK FOR TERMINATION CONDITIONS (trajectory specific)
    if !equilibrate_round && !do_adaptation()
        println("No more useful adaptations can be found.")
        break
    end

    # RUN OPTIMIZATION WITH CURRENT ANSATZ
    global converged = do_optimization()
    global equilibrate_round = false
end

!converged && println("Optimization did not converge. Run again to resume.")

JOB.save()
JOB.report()

name = "final"
JOB.archive(name)
JOB.plot(; label=name, trace=false)
JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
