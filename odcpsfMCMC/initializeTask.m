function task = initializeTask(task)
upSampleFactor = 4;

% initialize random number generator
task.randState = task.randSeed;    

% increase simulation region and setup simulation
if isfield(task,'dataGE')
    task.dataGE   = mypadarray(task.dataGE,size(task.dataGE),'post');
    task.dataSE   = mypadarray(task.dataSE,size(task.dataSE),'post');
    task.roiGE    = mypadarray(task.roiGE,size(task.roiGE),'post');
    task.roiSE    = mypadarray(task.roiSE,size(task.roiSE),'post');

    nVoxels1 = size(task.dataGE,1);
    nVoxels2 = size(task.dataGE,2);
else
    task.data   = mypadarray(task.data,size(task.data),'post');
    task.roi    = mypadarray(task.roi,size(task.roi),'post');

    nVoxels1 = size(task.data,1);
    nVoxels2 = size(task.data,2);
end

task.increaseFactor = 2;

vSize = 0.5;
task.sim = setupsim(nVoxels1,nVoxels2,vSize,upSampleFactor);

% initialize state
if isfield(task,'dataGE')
    noise0  = dcupsample(task.dataGE+task.dataSE,upSampleFactor);
else
    noise0  = dcupsample(task.data,upSampleFactor);
end
noise0  = 0.1 * noise0/std(noise0(:));
rho0    = 0.57;
delta0  = 0.24;
epsilon0= 0.67;
omega0  = (0.3 + 2)/2;
theta0  = 2.3;

if isfield(task,'dataGE')
    fwhmGE0 = 2;
    fwhmSE0 = 2;
else
    fwhm0   = 2;
end

if isfield(task,'dataGE')
    task.x = parametersToVectorSimultaneous(...
        noise0,rho0,delta0,epsilon0,omega0,theta0,fwhmGE0,fwhmSE0);
else
    task.x = parametersToVector(...
        noise0,rho0,delta0,epsilon0,omega0,theta0,fwhm0);
end

task.g = [];
task.E = [];
task.r = 0.651;

% setup function handle
if isfield(task,'dataGE')
    task.f = @(x) modelFSimultaneous(x,task.sim,task.dataGE,task.dataSE,...
        task.roiGE,task.roiSE,task.noiseLevelGE,task.noiseLevelSE,...
        task.betaGE,task.betaSE,...
        task.filterCutoffsGE,task.filterCutoffsSE,task.limitVector);
else
    task.f = @(x) modelF(x,task.sim,task.data,task.roi,...
        task.noiseLevel,task.beta,task.filterCutoffs,task.limitVector);
end

% setup MC parameters
task.nTau = 20;
task.epsilon = 0.005;

% initialize sampling
task.samples = zeros(2 * task.nsamples,length(task.x));
task.L = zeros(2 * task.nsamples,1);
task.zTrials = 0;
task.zSamples = 0;
end
