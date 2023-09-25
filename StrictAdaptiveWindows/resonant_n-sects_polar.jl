#= resonant_uniform

Resonant: Fix pulse frequencies on resonance.
n-sects: Split a single window into a rational division at each adaption.
Polar: Use complex pulses with polar parameterization.

=#

import StrictAdaptiveWindows as JOB
import StrictAdaptiveWindows: _!

import CtrlVQE

import LinearAlgebra: eigen, Hermitian, norm

using Serialization: serialize, deserialize
using NPZ: npzread
import SciPy

##########################################################################################
#= INITIALIZE SETUP VARIABLES =#

# LOAD THE MOLECULAR SYSTEM (INCLUDING H, REFERENCE STATE, EXACT GROUNDSTATE)
system = JOB.MolecularSystems.load_system(_!.setup.systemcode)

# PREP EVOLUTION ALGORITHM
evolution = CtrlVQE.Toggle(_!.setup.r)

# PREP DEVICE - THE PULSE IS A PLACEHOLDER UNTIL CURRENT STATE IS LOADED
protopulse = CtrlVQE.PolarComplexConstant(zero(JOB.Float), zero(JOB.Float))     # POLAR
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
        Ω = collect(1:2:L), φ = collect(2:2:L), ν = [],         # POLAR RESONANT
        s = [device.Ω̄[i].starttimes for i in 1:CtrlVQE.ndrives(device)],
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

cnt_f = isempty(_!.trace.f_calls) ? Ref(0) : Ref(last(_!.trace.f_calls))
lst_f = Ref(zero(JOB.Float))
function fun(energyfn)
    f = CtrlVQE.cost_function(energyfn);
    return x -> (
        value = f(x);
        lst_f[] = value;
        cnt_f[] += 1;
        value
    )
end

cnt_g = isempty(_!.trace.g_calls) ? Ref(0) : Ref(last(_!.trace.g_calls))
lst_g = Ref(zero(JOB.Float))
function jac(energyfn)
    g = CtrlVQE.grad_function(energyfn);
    return x -> (
        grads = g(x);
        lst_g[] = norm(grads, Inf);
        cnt_g[] += 1;
        grads
    )
end

iteration = isempty(_!.trace.iterations) ? Ref(0) : Ref(last(_!.trace.iterations))
function callback(x)
    iteration[] += 1

    # UPDATE TRACES
    push!(_!.trace.iterations, iteration[])
    push!(_!.trace.f_calls,    cnt_f[])
    push!(_!.trace.g_calls,    cnt_g[])

    push!(_!.trace.fn,         lst_f[])
    push!(_!.trace.gd,         lst_g[])

    # UPDATE DATA
    _!.state.x .= x
    JOB.save()
    if iteration[] % _!.opt.update == 0
        JOB.report()
    end

    # CHECK FOR SPECIAL TERMINATION CONDITIONS
    iteration[] > 10 && cnt_f[] > _!.setup.fnRATIO * iteration[] && (
        println("Linesearch excessively difficult");
        return true;
    )
    return false
end

function bounds()
    # PREPARE BOUND VECTORS
    L = length(_!.state.x)
    #= I thought that there was some documentation which said that
        L-BFGS-B needs boundless parameters to take "None" rather than "Inf".
    Alas, SciPy.jl (or maybe scipy, I dunno) seems to convert the bounds to float,
        so I *can't* pass `nothing`.
    Fortunately, using ±Inf seems to work anyways. =#
    μL = Vector{JOB.Float}(undef, L);
            μR = Vector{JOB.Float}(undef, L)
    μL[_!.state.Ω] .= -_!.setup.ΩMAX;                           # POLAR
            μR[_!.state.Ω] .= _!.setup.ΩMAX
    μL[_!.state.φ] .= -Inf;
            μR[_!.state.φ] .= Inf
    μL[_!.state.ν] .= xi[_!.state.ν] .- _!.setup.ΔMAX;
            μR[_!.state.ν] .= xi[_!.state.ν] .+ _!.setup.ΔMAX
    return collect(zip(μL, μR))
end

options = Dict(
    :iprint => 1,
    :ftol => _!.setup.f_tol,        # TOLERANCE IN FUNCTION VALUE
    :gtol => _!.setup.g_tol,        # TOLERANCE IN GRADIENT NORM
    :maxiter => _!.setup.maxiter,   # MAXIMUM NUMBER OF ITERATIONS
)

function do_optimization()
    optimization = SciPy.optimize.minimize(
        fun(energyfn),
        copy(_!.state.x);
        method="L-BFGS-B",
        jac=jac(energyfn),
        bounds=bounds(),
        options=options,
        callback=callback,
    )

    _!.state.x .= optimization["x"]
    return optimization["success"]
end

##########################################################################################
#= DEFINE ADAPTATION PROTOCOL =#

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

                # CALCULATE CHANGE IN GRADIENT
                g = CtrlVQE.Devices.gradient(candidate, τ̄, t̄, ϕ)
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

    # UPDATE THE DEVICE PULSES
    adapted = false
    # for i in 1:CtrlVQE.ndrives(device)    # TETRIS LOOP
    for i in argmax(ΔgBEST)                 # OPTIMAL "LOOP" (only one iteration)
        sBEST[i] == 0 && continue   # NO s GAVE SUFFICIENT IMPROVEMENT FOR THIS PULSE
        adapted = true              # *SOME* CHANGE WAS MADE IN THIS ADAPTATION

        device.Ω̄[i] = JOB.adapt_windows(device.Ω̄[i], sBEST[i])
    end
    !adapted && return false

    # UPDATE THE STATE
    x = CtrlVQE.Parameters.values(device)
    L = length(x)
    _!.state = JOB.StateVariables(
        x=x,
        Ω = collect(1:2:L), φ = collect(2:2:L), ν = [],         # POLAR RESONANT
        s = [device.Ω̄[i].starttimes for i in 1:CtrlVQE.ndrives(device)],
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
