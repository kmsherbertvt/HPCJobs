#= resonant_uniform

Resonant: Fix pulse frequencies on resonance.
Complex: Use complex pulses with inscribed amplitude bound.
Smooth: Use smooth bounds for amplitude and drive frequency.

=#

import BFGSUniformWindows as JOB
import BFGSUniformWindows: _!

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
protopulse = CtrlVQE.ComplexConstant(zero(JOB.Float), zero(JOB.Float))  # COMPLEX
device = CtrlVQE.Systematic(
    CtrlVQE.FixedFrequencyTransmonDevice,
    system.n,
    CtrlVQE.UniformWindowed(protopulse, _!.setup.T, _!.setup.W);
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
        JOB.archive(JOB.iterid(iteration[]))
        JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
    end

    # CHECK FOR SPECIAL TERMINATION CONDITIONS
    iteration[] > 10 && lossfn.f_counter[] > _!.opt.fnRATIO * iteration[] && (
        println("Linesearch excessively difficult");
        return true;
    )
    return false
end

options = Optim.Options(
    f_tol = _!.opt.f_tol,
    g_tol = _!.opt.g_tol,
    iterations = _!.opt.maxiter,
    callback = callback,
    store_trace = true,             # AFAIK, THIS IS THE ONLY WAY TO ACCESS Hk
    extended_trace = true,
)

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
#= RUN OPTIMIZATION =#

!_!.run && error("Setup Finished")

name = "init"
JOB.archive(name)
JOB.plot(; label=name, trace=false)

!do_optimization() && println("Optimization did not converge. Run again to resume.")

JOB.save()
JOB.report()

name = "final"
JOB.archive(name)
JOB.plot(; label=name, trace=false)
JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
