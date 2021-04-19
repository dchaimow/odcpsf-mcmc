# odcpsf-mcmc
Matlab code to fit a model of ODCs to fMRI data using MCMC as used in:

Chaimow, D., Yacoub, E., Uğurbil, K. & Shmuel, A. Spatial specificity of the functional MRI blood oxygenation response relative to neuronal activity. Neuroimage 164, 32–47 (2018).

`processMCMC.m` - script that demonstrates how to run MCMC on processed ODC fMRI data (makes some assumptions as to how fMRI results are saved, see [loadDataForFitting.m](+runMCMC/private/loadDataForFitting.m))

`odcpsfMCMC` - contains the main code of the model and of the MCMC algorithm

`+runMCMC` - contains code involved in starting and running MCMC jobs locally or on an HPC cluster

`+tests` - contains tests

`+tools` - contains some helper functions




