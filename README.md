# odcpsf-mcmc
Matlab code to fit a model of ODCs to fMRI data using MCMC as used in:

Chaimow, D., Yacoub, E., Uğurbil, K. & Shmuel, A. Spatial specificity of the functional MRI blood oxygenation response relative to neuronal activity. Neuroimage 164, 32–47 (2018).

`processMCMC.m` - script that demonstrates how to run MCMC on processed ODC fMRI data (makes some assumptions as to how fMRI results are saved, see [loadDataForFitting.m](+runMCMC/private/loadDataForFitting.m))

`odcpsfMCMC` - the main code of the model and of the MCMC algorithm

`+runMCMC` - code involved in starting and running MCMC jobs locally or on an HPC cluster

`+tests` - tests

`+tools` - some helper functions

## General description

The generative model of fMRI of ODCs is implemented in [computeODCModel.m](odcpsfMCMC/computeODCModel.m) (and [computeODCModelSimultaneous.m](odcpsfMCMC/computeODCModelSimultaneous.m)). This function simulates ODC fMRI maps given a set of parameters and also computes the energy and the gradient given real fMRI data.


