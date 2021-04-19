# odcpsf-mcmc
Matlab code to fit a model of ODCs to fMRI data using MCMC as used in:

Chaimow, D., Yacoub, E., Uğurbil, K. & Shmuel, A. Spatial specificity of the functional MRI blood oxygenation response relative to neuronal activity. Neuroimage 164, 32–47 (2018).

`processMCMC.m` - script that demonstrates how to run MCMC on processed ODC fMRI data (makes some assumptions as to how fMRI results are saved, see [loadDataForFitting.m](+runMCMC/private/loadDataForFitting.m)) and visualizes results.

`odcpsfMCMC` - the main code of the model and of the MCMC algorithm

`+runMCMC` - code involved in starting and running MCMC jobs locally or on an HPC cluster

`+tests` - tests

`+tools` - some helper functions

## General description

The generative model of fMRI of ODCs is implemented in [computeODCModel.m](odcpsfMCMC/computeODCModel.m) (and [computeODCModelSimultaneous.m](odcpsfMCMC/computeODCModelSimultaneous.m)). This function simulates ODC fMRI maps given a set of parameters and also computes the energy and the gradient given real fMRI data.

The Hamiltonian Monte Carlo algorithm is implemented independent of the model in [dcHMCRunOnce.m](odcpsfMCMC/dcHMCRunOnce.m).

Generally, in order to run MCMC model fitting on real data:

1. Create an array of `task` structures from processed fMRI data and desired options using [createMCFittingTasks.m](+runMCMC/createMCFittingTasks.m), which in turn uses [loadDataForFitting.m](+runMCMC/private/loadDataForFitting.m) to load processed fMRI data (also does some outlier processig). Alternatively `task` could be set manually by setting:
   - `task.data` - differential response map
   - `task.roi` - boolean roi mask
   - `task.beta` - response level, estimated as twice the median of responses averaged between both conditions
   - `task.noiseLevel` - noise level as RMS of standard error of differential responses
   - `task.filterCutoffs` - bandpass filter cutoffs, if used in processing of fMRI data
   - `task.noBVMaskFlag` - true if blood vessels voxels were not masked out
   - `task.limitVector` - cell array of parameter estimation limits for rho and omega
   - `task.removeMean` - true if mean was removed from differential maps
   - `task.nsamples` - number of samples to obtain  (After reaching this number of samples, every 2nd sample is removed = thin level 2, sampling continues recording every 2nd sample until nsamples is reached again, after which only every 4th sample is recorded = thin level 4. This process repeats until a set time or thin level is reached.)
   - `task.thin` - thin level to start with, usually 1
   - `task.randSeed` - seed for random number generator
   - `task.sesname` - descriptor for processed data set
   - `task.taskID` - task ID
2. Create a job on a HPC cluster or locally using [runMCFittingTasks.m](+runMCMC/runMCFittingTasks.m) and let it run. Local tasks can also be run directly using runMCMC.runSingleMCMCTaskLocal(tasks). There are a number of script that are involved in compiling and uploading the code to the HPC cluster, to restart jobs and to downloa results.  
3. Running MCMC on the model basically means that (e.g. see [runSingleMCMCTaskLocal](+runMCMC/runSingleMCMCTaskLocal.m)):
   1. The task structure is initialized by passing it to  [initializeTask.m](odcpsfMCMC/dcHMCRunOnce.m) in order to prepare it for MCMC sampling
   2. The HMC sampling algorithm is called iteratively with a function handle that has been generated in i. by wrapping [computeODCModel.m](odcpsfMCMC/computeODCModel.m) using [modelF.m](odcpsfMCMC/modelF.m) with the current set of parameters in a parameter vector and the current energy and gradient value. It will then sample a new parameter vector and also return the new energy and gradient 
