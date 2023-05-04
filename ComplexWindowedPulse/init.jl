#= Read in parameters and initialize a job directory.

# Usage
```julia init.jl matrix T W r m ```

# Parameters
- matrix: filename (without path or extension) of matrix
- T: Pulse duration in ns
- W: Number of windows in the square pulse.
- r: Number of Trotter steps.
- m: Size of Hilbert space for each qubit.

For convenience,
    omitting m defaults to 2,
    omitting r defaults to 20T, and
    omitting W defaults to 5T.
You can't omit earlier arguments without omitting later arguments.

=#

shift!(array, default=nothing) = length(array) > 0 ? popfirst!(array) : default

# CONSTRUCT jobdir AND ENSURE IT IS AVAILABLE
jobdir = "CWP_"*join(ARGS, "_")
isdir(jobdir) && error("Some data already exists! Delete or rename directory '$jobdir'")

# PARSE ARGUMENTS AS SETTINGS
import ComplexWindowedPulse as JOB

matrix = shift!(ARGS)
T = parse(Float64, shift!(ARGS))
W = shift!(ARGS, nothing)
r = shift!(ARGS, nothing)
m = shift!(ARGS, nothing)

JOB.modify_settings(
    matrixfile="HPCJobs/matrix/$matrix.npy",
    T=T,
    W= !isnothing(W) ? parse(Int, W) : round(Int,10T),
    r= !isnothing(r) ? parse(Int, r) : round(Int,20T),
    m= !isnothing(m) ? parse(Int, m) : 2,
)
mkpath(jobdir)
JOB.save_settings(jobdir)

println(jobdir)     # EXPORTS `jobdir` TO THE CALLING BASH SCRIPT