#= resonant_uniform

Resonant: Fix pulse frequencies on resonance.
Complex: Use complex pulses with inscribed amplitude bound.
Smooth: Use smooth bounds for amplitude and drive frequency.

=#

import StrictUniformWindows as JOB
import StrictUniformWindows: _!

import CtrlVQE

import LinearAlgebra: eigen, Hermitian, norm, diagm

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
        JOB.archive(JOB.iterid(iteration[]))
        JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
    end

    # CHECK FOR SPECIAL TERMINATION CONDITIONS
    iteration[] > 10 && cnt_f[] > _!.opt.fnRATIO * iteration[] && (
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
    μL[_!.state.Ω] .= -_!.setup.ΩMAX / √2;              # COMPLEX (INSCRIBED)
            μR[_!.state.Ω] .= _!.setup.ΩMAX / √2
    μL[_!.state.φ] .= -Inf;
            μR[_!.state.φ] .= Inf
    μL[_!.state.ν] .= xi[_!.state.ν] .- _!.setup.ΔMAX;
            μR[_!.state.ν] .= xi[_!.state.ν] .+ _!.setup.ΔMAX
    return collect(zip(μL, μR))
end

options = Dict(
    :iprint => 1,
    :ftol => _!.opt.f_tol,          # TOLERANCE IN FUNCTION VALUE
    :gtol => _!.opt.g_tol,          # TOLERANCE IN GRADIENT NORM
    :maxiter => _!.opt.maxiter,     # MAXIMUM NUMBER OF ITERATIONS
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
#= RUN OPTIMIZATION =#

!_!.run && error("Setup Finished")

name = "init"
JOB.archive(name)
JOB.plot(; label=name, trajectory=false, trace=false)

!do_optimization() && println("Optimization did not converge. Run again to resume.")

JOB.save()
JOB.report()

name = "final"
JOB.archive(name)
JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
JOB.plot(; label=name, trace=false) # INCLUDE TRAJECTORY
