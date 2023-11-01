# HPCJobs
Facilitates transfering files between desktop and remote HPC resources

Principles:
- Low-resource job initialization that puts everything a job needs to run in a job folder.
- Simple interface for starting a job, and almost-as-easy interface to do so with slurm.
- Jobs are NOT saved to git. Transfer data files with scp.
- Packages provide tools to _inspect_ job data _without_ any expensive calculations.

Preferrred Workflow:
- Package defines state variables in an object, which can be archived throughout a job.
- Package defines an easy `load` function to load any archived state as the "active" one.
- The package defines a global variable `!_` to hold the active state.
- Plotting and other analysis methods act directly on the active state.

I expect this business with the `!_` variable will look very odd to any visitors,
    but I've found it works extraordinarily well for me.