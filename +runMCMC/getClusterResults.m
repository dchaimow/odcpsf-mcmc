function results = getClusterResults(jobs,descriptionString)
nJobs = length(jobs);
resultsPerJob = cell(1,nJobs);
% load results per job
for zJob=1:nJobs
    if isa(jobs,'parallel.job.CJSCommunicatingJob')
        fname = jobs(zJob).fetchOutputs{1};
        if ~strcmp(jobs(zJob).Parent.Profile,'local')
            system(['rsync -avz --remove-source-files ' ...
                'dchaimow@guillimin.clumeq.ca:' fname '.mat .']);
        end
        load(fname);
        resultsPerJob{zJob} = results;
    else
        system(['rsync -avz --remove-source-files ' ...
            'dchaimow@guillimin-p2.hpc.mcgill.ca:' ...
            fullfile('/home/dchaimow/jobs/',jobs{zJob}) ' .']);
        load(fullfile(jobs{zJob},'results.mat'));
        resultsPerJob{zJob} = results;
    end
end
% merge all results
results = cat(2,resultsPerJob{:});

% sort all results
[~, sortOrder] = sort([results(:).taskID],'ascend');
results = results(sortOrder);

% save results
if ~exist('descriptionString','var')
    descriptionString = '';
else
    descriptionString = [descriptionString '_'];
end
fname = fullfile('~/Data/psfodc/results',...
    ['odcpsf_MCMC_results_' ...
    descriptionString datestr(now,30)]);
save(fname,'results','-v7.3');

% remove local files
for zJob=1:nJobs
    if isa(jobs,'parallel.job.CJSCommunicatingJob')
        fname = jobs(zJob).fetchOutputs{1};
        delete([fname '.mat']);
    else
        rmdir(jobs{zJob},'s');
    end
end
