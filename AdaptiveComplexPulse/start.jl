shift!(array, default=nothing) = length(array) > 0 ? popfirst!(array) : default

jobdir = shift!(ARGS)
!isdir(jobdir) && error("$jobdir/ must exist. Consider running 'julia init.jl ...'")
mode = shift!(ARGS, nothing)

# SELECT DEFAULT BASED ON STATE OF DIRECTORY (use as much data as is available)
if isnothing(mode)
    mode = (
        isfile("$jobdir/iterations") ?
            "--from-where-we-left-off" :
        isfile("$jobdir/x") ?
            "--from-x" :
            "--from-init"
    )
end

import AdaptiveComplexPulse as JOB
JOB.load_settings(jobdir)
if mode == "--from-init"
    JOB.load_settings(jobdir)
    JOB.initialize()
elseif mode == "--from-x"
    JOB.load(jobdir, resume=false)
else
    JOB.load(jobdir, resume=true)
end

JOB.run_optimization(jobdir)
JOB.update(jobdir)

while JOB.run_adaptation(jobdir)
    JOB.run_optimization(jobdir)
    JOB.update(jobdir)
end