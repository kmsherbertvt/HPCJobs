import AdaptiveWindows as JOB
import AdaptiveWindows: _!, Float

import CtrlVQE

import Random: seed!
import LinearAlgebra: eigen, Hermitian

import NPZ: npzread
import Optim, LineSearches

##########################################################################################
#= INITIALIZE SETUP VARIABLES =#

H = npzread("$(ENV["HOME"])/HPCJobs/matrix/$(_!.setup.matrix).npy")
n = CtrlVQE.QubitOperators.nqubits(H)
ψ_REF = CtrlVQE.QubitOperators.reference(H)
REF = real(ψ_REF' * H * ψ_REF)
Λ, U = eigen(Hermitian(H))

pulse = CtrlVQE.Signals.WindowedSignal(
    [CtrlVQE.Signals.ComplexConstant(zero(Float), zero(Float))],
    [0.0],
)
systematicdevice = CtrlVQE.SystematicTransmonDevice(_!.setup.m, n, pulse)
device = CtrlVQE.Devices.FixedFrequencyTransmonDevice(
    systematicdevice.ω̄,
    systematicdevice.δ̄,
    systematicdevice.ḡ,
    systematicdevice.quples,
    systematicdevice.q̄,
    systematicdevice.ν̄,
    systematicdevice.Ω̄,
    systematicdevice.m,
)
algorithm = CtrlVQE.Rotate(_!.setup.r)

τ, _, t = CtrlVQE.trapezoidaltimegrid(_!.setup.T, _!.setup.r)

seed!(_!.setup.seed)
xi = CtrlVQE.Parameters.values(device)

L = length(xi)
# Ω = 1:L-n; φ = []; ν = 1+L-n:L              # Cartesian with Frequencies
Ω = 1:L; φ = []; ν = []                     # Cartesian without Frequencies

xi[Ω] .+= _!.setup.init_Ω .* (2 .* rand(length(Ω)) .- 1)
xi[φ] .+= _!.setup.init_φ .* (2 .* rand(length(φ)) .- 1)
xi[ν] .+= _!.setup.init_Δ .* (2 .* rand(length(ν)) .- 1)

O0 = CtrlVQE.QubitOperators.project(H, device)
ψ0 = CtrlVQE.QubitOperators.project(ψ_REF, device)

fn_energy, gd_energy = CtrlVQE.ProjectedEnergy.functions(
    O0, ψ0, _!.setup.T, device, _!.setup.r;
    frame=CtrlVQE.STATIC,
)

fn_norm, gd_norm = CtrlVQE.Normalization.functions(ψ0, _!.setup.T, device, _!.setup.r)



##########################################################################################
#= INITIALIZE OPTIMIZATION VARIABLES =#

function callback(state)
    # EXTRACT STATE OF *CURRENT* OPTIMIZATION
    iteration = state.iteration
    f_calls = fn.counter[]
    g_calls = gd.counter[]

    # UPDATE TRACES
    push!(_!.trace.iterations, iteration + offset_iterations)
    push!(_!.trace.f_calls,    f_calls   + offset_f_calls)
    push!(_!.trace.g_calls,    g_calls   + offset_g_calls)

    push!(_!.trace.fn,         state.value)
    push!(_!.trace.fn_energy,  fn.values[1])
    push!(_!.trace.gd,         state.g_norm)
    push!(_!.trace.gd_energy,  gd.norms[1])
    push!(_!.trace.gd_penalty, gd.norms[2])

    # UPDATE DATA
    _!.state.x .= CtrlVQE.Parameters.values(device)
    JOB.save()
    if iteration % _!.opt.update == 0
        JOB.report()
    end

    # CHECK FOR SPECIAL TERMINATION CONDITIONS
    iteration > 10 && f_calls > _!.setup.fnRATIO * iteration && (
        println("Linesearch excessively difficult");
        return true;
    )
    return false
end

if length(_!.trace.iterations) > 0
    offset_iterations = last(_!.trace.iterations)
    offset_f_calls = last(_!.trace.f_calls)
    offset_g_calls = last(_!.trace.g_calls)
else
    offset_iterations = 0
    offset_f_calls = 0
    offset_g_calls = 0
end

linesearch = LineSearches.MoreThuente()

optimizer = Optim.LBFGS(
    linesearch=linesearch,
)

options = Optim.Options(
    show_trace = true,
    show_every = 1,
    f_tol = _!.setup.f_tol,
    g_tol = _!.setup.g_tol,
    iterations = _!.setup.maxiter,
    callback = callback,
)


##########################################################################################
#= INITIALIZE STATE VARIABLES (if necessary) =#
if isnothing(_!.state)
    _!.state = JOB.StateVariables(
        x = copy(xi),
        Ω = collect(Ω),
        φ = collect(φ),
        ν = collect(ν),
        s = [[0.0] for i in 1:CtrlVQE.ndrives(device)],
    )
else
    cursor = 1
    for i in 1:CtrlVQE.ndrives(device)
        windows = CtrlVQE.Signals.ComplexConstant{Float}[]
        for (k, s) in enumerate(_!.state.s[i])
            push!(windows, CtrlVQE.Signals.ComplexConstant(
                _!.state.x[cursor],
                _!.state.x[cursor+1],
            ))
            global cursor += 2
        end
        device.Ω̄[i] = CtrlVQE.Signals.WindowedSignal(windows, _!.state.s[i])
    end
end
CtrlVQE.Parameters.bind(device, _!.state.x)

##########################################################################################
#= ASSIGN WORK VARIABLES =#

_!.work = JOB.WorkVariables(
    REF = REF,
    FCI = Λ[1],
    FES = Λ[2],

    ψ_FCI = U[:,1],

    device = device,
    fn_energy = fn_energy,
    gd_energy = gd_energy,
    fn_norm = fn_norm,
)


##########################################################################################
#= DEFINE OPTIMIZATION AND ADAPTATION PROTOCOLS =#

function do_optimization()
    L = length(_!.state.x)

    λ  = zeros(L);  λ[_!.state.Ω] .=    _!.setup.λΩ
    μR = zeros(L); μR[_!.state.Ω] .= +_!.setup.ΩMAX ./ √2
    μL = zeros(L); μL[_!.state.Ω] .= -_!.setup.ΩMAX ./ √2
    σ  = zeros(L);  σ[_!.state.Ω] .=    _!.setup.σΩ
    fn_penalty, gd_penalty = CtrlVQE.SmoothBounds.functions(λ, μR, μL, σ)

    global fn = CtrlVQE.CompositeCostFunction(_!.work.fn_energy, fn_penalty)
    global gd = CtrlVQE.CompositeGradientFunction(_!.work.gd_energy, gd_penalty)

    optimization = Optim.optimize(fn, gd, _!.state.x, optimizer, options)
    println(optimization)

    global offset_iterations = last(_!.trace.iterations)
    global offset_f_calls = last(_!.trace.f_calls)
    global offset_g_calls = last(_!.trace.g_calls)

    return Optim.converged(optimization)
end

function do_adaptation()
    _, τ̄, t̄ = CtrlVQE.trapezoidaltimegrid(_!.setup.T, _!.setup.r)

    # CALCULATE GRADIENT SIGNAL
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
    ϕ̄ = _!.work.gd_energy.ϕ̄

    # CALCULATE OPTIMAL SWITCHING TIMES FOR EACH PULSE
    ΔgBEST = zeros(CtrlVQE.ndrives(_!.work.device))  # TRACK LARGEST Δg FOR EACH PULSE
    sBEST  = zeros(CtrlVQE.ndrives(_!.work.device))  # TRACK s GIVING LARGEST Δg " "

    maxnsect = floor(Int, _!.setup.T / _!.setup.ΔsMIN)
    for nsect in 2:maxnsect
        i_tosplit = findall(s -> s == 0, sBEST) # WHICH DRIVES STILL NEED A SPLIT?
        isempty(i_tosplit) && break             # STOP ONCE ALL DRIVES HAVE A SPLIT

        for j in 1:CtrlVQE.ngrades(_!.work.device)
            i = (j-1)÷2 + 1                     # IDENTIFY DRIVE
            !(i in i_tosplit) && continue       # SKIP DRIVE SPLITTABLE WITH SMALLER nsect

            pulse = _!.work.device.Ω̄[i]
            s̄ = JOB.nsect_windows(nsect, pulse.starttimes)
            JOB.filter_splits!(s̄, pulse.starttimes)

            for s in s̄
                # CONSTRUCT TRIAL DEVICE
                candidate = deepcopy(device)
                candidate.Ω̄[i] = JOB.adapt_windows(device.Ω̄[i], s)

                # CALCULATE CHANGE IN GRADIENT
                g = CtrlVQE.Devices.gradient(candidate, τ̄, t̄, ϕ̄)
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
    # for i in 1:CtrlVQE.ndrives(_!.work.device)    # TETRIS LOOP
    for i in argmax(ΔgBEST)
        sBEST[i] == 0 && continue   # NO s GAVE SUFFICIENT IMPROVEMENT FOR THIS PULSE
        adapted = true              # *SOME* CHANGE WAS MADE IN THIS ADAPTATION

        _!.work.device.Ω̄[i] = JOB.adapt_windows(device.Ω̄[i], sBEST[i])
    end
    !adapted && return false

    # UPDATE THE STATE
    s = [_!.work.device.Ω̄[i].starttimes for i in 1:CtrlVQE.ndrives(_!.work.device)]
    x = CtrlVQE.Parameters.values(_!.work.device)
    L = length(x)
    # Ω = 1:L-n; φ = []; ν = 1+L-n:L              # Cartesian with Frequencies
    Ω = 1:L; φ = []; ν = []                     # Cartesian without Frequencies
    _!.state = JOB.StateVariables(x=x, Ω=Ω, φ=φ, ν=ν, s=s)

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

converged = do_optimization()
if !converged
    println("Optimization did not converge. Run again to resume.")
else
    JOB.save()
    JOB.report()

    name = "optimized_$(JOB.adaptid(length(_!.trace.adaptations)))"
    JOB.archive(name)
    JOB.plot(; label=name, trace=false)
    JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
end

while converged
    if length(_!.trace.adaptations) ≥ _!.opt.maxadapt
        println("Maximum adaptations reached.")
        break
    end

    if !do_adaptation()
        println("No more useful adaptations can be found.")
        break
    end

    if !do_optimization()
        println("Optimization did not converge. Run again to resume.")
        break
    else
        JOB.save()
        JOB.report()

        name = "optimized_$(JOB.adaptid(length(_!.trace.adaptations)))"
        JOB.archive(name)
        JOB.plot(; label=name, trace=false)
        JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
    end
end

JOB.save()
JOB.report()

name = "final"
JOB.archive(name)
JOB.plot(; label=name, trace=false)
JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
