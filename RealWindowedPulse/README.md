# Workflow

To start a new job, `cd` into this project diretory and then
```
./run {matrix} {T}
```

To see the fruits of your labor,
```
./inspect {jobdir}
```

To resume an optimization with a lower tolerance,
```
julia
: jobdir = "..." # or whatever
: ]activate .
: import RealWindowedPulse as JOB
: JOB.load({jobdir})
: JOB.modify_settings(
:     g_tol=1e-8, # or whatever
: )
: JOB.save(jobdir)
: exit()
./run {jobdir} --from-x
```

