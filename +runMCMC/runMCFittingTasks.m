function jobs = runMCFittingTasks(tasks,nLocalCores)
% jobs = runMCFittingTasks(tasks,nLocalCores)
%
% runs locally on nLocalCores cores
% of nLocalCores not supplied, runs on Guillimin
% returns job structure
% retrieve data after job is done with:
% results = getClusterResults(jobs)
addpath('../odcpsfMCMC');

nTasks = length(tasks);
% randomize task order
tasks = tasks(randperm(nTasks));

% set cluster and distribute tasks onto jobs
if exist('nLocalCores','var')
    % on local cluster start one job on a nLocalClores-1 matlabpool
    % resulting in nLocalCores being used
    myCluster = parcluster('local');
    tasksPerJob = {tasks};
    nTasksPerJob = nTasks;
    mpoolSize = nLocalCores-1;
    nJobs = 1;
    sampleTime = 1714/1000;
else
    % on hpc cluster each job is limited to one node with 12 cores
    % distribute tasks across multiple jobs so that each task
    % uses one core; 1 extra core is needed per job
    myCluster = parcluster('Guillimin');
    maxCores = 64;
    maxTasksPerJob = maxCores - 1;
    nJobs = ceil(nTasks/maxTasksPerJob);
    [tasksPerJob, nJobs] = pctdemo_helper_split_vector(tasks, nJobs);    
    nTasksPerJob = cellfun(@length,tasksPerJob); 
    mpoolSize = nTasksPerJob;
    sampleTime = 1690/1000;
end

% estimate running time
nSamplesPerTask = (sum([tasks(:).nsamples])/nTasks) ...
    *  (isfield(tasks(:),'dataGE')+1);
nThinPerTask = sum([tasks(:).thin])/nTasks;
singleTaskTime = nThinPerTask * nSamplesPerTask * sampleTime;
estimatedRunningTime = (singleTaskTime .* nTasksPerJob)./mpoolSize;

% sumbit jobs
nMaxProcessorsPerNode = 16;
jobs = parallel.job.CJSCommunicatingJob.empty(0,nJobs);
for zJob = 1:nJobs
    nNodes = ceil((mpoolSize(zJob) + 1)/nMaxProcessorsPerNode);
    nProcessorsPerNode = ceil((mpoolSize(zJob) + 1)/nNodes);
    writeMoabFile(nNodes,nProcessorsPerNode,...
        max(estimatedRunningTime(zJob)*2,10*60));
    disp(['job ' num2str(zJob) '/' num2str(nJobs) ...
        ' , estimated running time: ' ...
        datestr(estimatedRunningTime(zJob)/(24*60*60),'DD:HH:MM:SS')]);
    jobs(zJob) = batch(myCluster,@runFitting,1,tasksPerJob(zJob),...
        'pool', mpoolSize(zJob),'CurrentDirectory', '.');
end
end

function writeMoabFile(nodes,ppn,timeInSec)
timeStr = datestr(timeInSec/86400,'DD:HH:MM:SS');
str0 = 'function [nodes,ppn,time] = moab()\n';
str1 = ['nodes = ' num2str(nodes) ';\n'];
str2 = ['ppn = ' num2str(ppn) ';\n'];
str3 = ['time = ''' timeStr ''';\n'];
str4 = 'end';

str = [str0 str1 str2 str3 str4];

fid = fopen('moab.m','w');
fprintf(fid, str);
fclose(fid);
end

function fname = runFitting(tasks)
nTasks = length(tasks);
parfor zTask=1:nTasks
    results(zTask) = runMCMC.runSingleMCMCTaskLocal(tasks(zTask)); %#ok<NASGU,PFOUS>
end
if isfield(tasks(1),'jobString')
    jobString = [tasks(1).jobString '_'];
else
    jobString = '';
end
fname = ['odcpsf_MCMC_results_' jobString num2str(tasks(1).taskID) ...
    '_' datestr(now,30)];
save(fname,'results','-v7.3');
end
