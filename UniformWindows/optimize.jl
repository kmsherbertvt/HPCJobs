import UniformWindows as JOB
import UniformWindows: _!, Float

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

starttimes = range(zero(_!.setup.T), _!.setup.T, _!.setup.W+1)[1:end-1]
pulse = CtrlVQE.Signals.WindowedSignal(
    [CtrlVQE.Signals.ComplexConstant(zero(Float), zero(Float)) for t in starttimes],
    starttimes,
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
Ω = 1:L; φ = []; ν = []                   # Cartesian without Frequencies

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

λ  = zeros(L);  λ[Ω] .=    _!.setup.λΩ
μR = zeros(L); μR[Ω] .= +_!.setup.ΩMAX ./ √2
μL = zeros(L); μL[Ω] .= -_!.setup.ΩMAX ./ √2
σ  = zeros(L);  σ[Ω] .=    _!.setup.σΩ
fn_penalty, gd_penalty = CtrlVQE.SmoothBounds.functions(λ, μR, μL, σ)

fn = CtrlVQE.CompositeCostFunction(fn_energy, fn_penalty)
gd = CtrlVQE.CompositeGradientFunction(gd_energy, gd_penalty)

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
        JOB.archive(JOB.iterid(last(_!.trace.iterations)))
    end

    # CHECK FOR SPECIAL TERMINATION CONDITIONS
    iteration > 10 && f_calls > _!.opt.fnRATIO * iteration && (
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
    f_tol = _!.opt.f_tol,
    g_tol = _!.opt.g_tol,
    iterations = _!.opt.maxiter,
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
    )
end

##########################################################################################
#= ASSIGN WORK VARIABLES =#

Λ, U = eigen(Hermitian(H))
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
#= RUN OPTIMIZATION =#

!_!.run && error("Setup Finished")

optimization = Optim.optimize(
    x -> (print("\tf:"); @time fn(x)),
    (G,x) -> (print("\tg:"); @time gd(G,x)),
    _!.state.x,
    optimizer,
    options,
)

JOB.save()
JOB.report()

name = "final"
JOB.archive(name)
JOB.plot(; label=name, trace=false)
JOB.plot(; label="", trajectory=false, pulses=false, trace=true)
