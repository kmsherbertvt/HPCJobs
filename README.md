# HPCJobs
Facilitates transfering files between desktop and remote HPC resources


The current job writing is nice in a lot of ways, but not nice in others.

Good:
- Simple command to just run with different experimental parameters.
- Easy post-processing, and capability to restart, and modify parameters.

Bad:
- Very unreadable. Hard to find parts of code that you want to stare at.
- Jobs don't remember non-parametric modifications, eg. line-search.

Solution to the Bad is to write a proper *script*, rather than this mumbo-jumbo of methods.
We can still use a module setup to do standard things like loading and plotting,
    and defining the parameter framework.
    But the actual actions should take place in a script format probably.
The script itself gets saved within each job alongside the experimental parameters.

We can use a mutable struct instead of the module to hold parameters and variables.
On init we create parameters, and blank variables, and save script.
On load we load parameters, and state variables, and trace variables.
On start we run the saved script using active parameters and variables.

The script can't take mutable struct as arg.
It will need to refer to a global variable. Something like `_!`.
`_!` lives in the job module, alongside struct definition and load/plot logic.
Script starts with the line `import JOB: vars as _!`