# `odcpsfMCMC` - the main code of the model and of the MCMC algorithm

## Model implementation

### `computeODCModel.m`
Generative ODC model.
Simulates an ODC (and MRI)  pattern given a set of parameters. Can also return the value of the energy function given data as well as its gradient. Also expects a simulation structure sim, that is generated using setupsim.m and specifies the simulation grid.

### `computeODCModelSimultaneous.m`
Same as above for a pair of MRI images resulting from different BOLD PSF widths and response amplitudes.


### Helper functions
#### `makeFFT2Symmetric.m`
Returns the (conjugate) symemtric part of a matrix. Is used here to correct small asymmetries in spatial frequency representations either due to numerical reasons or the inherent slight asymmetry of the FFT.

#### `fft2DoubleFlip.m`
Used in  `makeFFT2Symmetric.m`.

#### `fft2DoubleFlipConj.m`
Used in  `makeFFT2Symmetric.m`.

## MCMC implementation

### `dcHMCRunOnce.m`
Implementation of Hamiltonian (Hybrid)  Monte Carlo. Expects a function handle to a function [E,g] = f(x), with x being a vector of parameters, and E and g being the energy and the gradient of the model.

## Functions that run MCMC using ODC Model (self contained code that would run on HPC cluster)

# make.sh
Compilation script, compiles code before running on HPC cluster.

# odcpsfmcmc.m
Main function started on HPC cluster, calls `runSingleMCMCTask.m`. 

# odcpsfmcmcMultiThreaded.m
Alternative version using multithreaded compilation.

# restart_jobs.sh
Shell script to restart aborted jobs in HPC cluster.

### `runSingleMCMCTask.m`
Runs MCMC task by first running `initializeTask.m`, then looping over `dcHMCRunOnce.m`, managing, saving and visualizing results in between.

### `initializeTask.m`
takes a task structure and prepares it for running MCMC on it by
- increasing simulation space to twice the data size (in order to avoid periodicity artefacts)
- setting up simulation grid using setupsim
- initializing starting parameter state
- setting MCMC slgorithm parameters
- setting up energy function handle to pass to MCMC algorithm
- initializing sampling variables and counters

### `modelF.m`
wraps computeODCModel.m to use it as an energy function as required by dcHMCRunOnce.m

### `modelFSimultaneous.m`
wraps computeODCModelSimultaneous.m to use it as an energy function as required by dcHMCRunOnce.m


### Helper functions
`mypadarray.m`
`parametersToVector.m`
`parametersToVectorSimultaneous.m`
`setupsim.m`
`vectorToParameters.m`
`vectorToParametersSimultaneous.m`
`dcupsample.m`

## Visualization functions

### `visualizeMCDependencies.m`

### `visualizeMCDependenciesSimultaneous.m`

### `visualizeMCGoF.m`

### `visualizeMCGoFSimultaneous.m`

### `visualizeMCHist.m`

### `visualizeMCHistSimultaneous.m`

### `visualizeMCSamples.m`

### `visualizeMCSamplesSimultaneous.m`

### `visualizeResults.m`

### `displayResults.m`

### `displayResultsForTask.m`

### Helper functions

#### `bwrLinCMap.mat`
Colormap used in visualizing ODC maps.

#### `dcprintfig.m`
Helper function to exports a figure.

#### `imagebw.m`
Display black and white map.
