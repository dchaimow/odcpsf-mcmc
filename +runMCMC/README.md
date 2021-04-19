# `+runMCMC` - code involved in starting and running MCMC jobs locally or on an HPC cluster

## `createMCFittingTasks.m`
Creates task structures for MCMC sampling using processed fMRI data set and specified options.

## `createMCFittingTasksSimultaneous.m`
Same as above for simultaneous fitting.

## `getClusterResults.m`
Download results from submitted jobs.

## `runMCFittingTasks.m`
Submits tasks to run locally using matlab parallel computing toolbox.

## `runSingleMCMCTaskLocal.m`
Run a single task locally (no job submission, MCMC loop take place directly within this script).

## `runTasks.m`
Runs tasks by submitting them to HPC cluster (uses submitTasksForGuilliminCompiled)

## `submitTasksForGuilliminCompiled.m`
Submits tasks to be run on Compute Canada McGill HPC cluster Guillimin (out of service now).

## `uploadAndCompile.m`
Uploads code to HPC cluster and compiles it.
