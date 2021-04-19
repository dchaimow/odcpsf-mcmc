%% first make sure all data is processed with the right parameters
%addpath ../../..;
%processSubject('s1');
%processSubject('s2');

addpath('odcpsfMCMC');

jobDescription = 'test';

%% run MCMC
%upload and compile on HPC cluster
%runMCMC.uploadAndCompile();

nSamples = 100; % 100 for testing, otherwise use 1000
nThin = 1;
rhoLimits = [0.37 0.77];
omegaLimits = [0.3 2];
noBVMask = true;
removeMean = false;

% single
jobNameSingle = ...
    [jobDescription '_single_' datestr(now,'yyyy_mm_dd')];
tasks.(jobNameSingle) = ...
    runMCMC.createMCFittingTasks(...
    {'s1_ge','s1_se'},noBVMask,removeMean,{rhoLimits,omegaLimits},nSamples,nThin);

% simultaneous
jobNameSimultaneous = ...
    [jobDescription '_simultaneous_' datestr(now,'yyyy_mm_dd')];
tasks.(jobNameSimultaneous) = ...
    runMCMC.createMCFittingTasksSimultaneous(...
    {'s1_ge','s1_se'},noBVMask,removeMean,{rhoLimits,omegaLimits},nSamples,nThin);


% run locally
nCores = 4;
jobSingle = ...
    runMCMC.runMCFittingTasks(tasks.(jobNameSingle),nCores);
jobSimultaneous = ...
    runMCMC.runMCFittingTasks(tasks.(jobNameSimultaneous),nCores);

% run on server
% tasks.(jobNameSingle) = ...
%     runMCMC.runTasks(tasks.(jobNameSingle),jobNameSingle);
% tasks.(jobNameSimultaneous) = ...
%     runMCMC.runTasks(tasks.(jobNameSimultaneous),jobNameSimultaneous);

% visualize results
mkdir figures;
resultsFilenames = dir('odcpsf_MCMC_results_1*.mat');
load(resultsFilenames(1).name,'results');
visualizeResults(results,'figures');
load(resultsFilenames(2).name,'results');
visualizeResults(results,'figures');