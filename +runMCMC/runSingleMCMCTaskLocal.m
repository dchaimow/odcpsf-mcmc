function result = runSingleMCMCTaskLocal(task)
addpath('odcpsfMCMC');
if ~isfield(task,'x')
    task = initializeTask(task);
end
rng(task.randState)     % set random number generator
while task.zSamples < task.nsamples
    % update state
    [task.x,task.g,task.E,task.r,task.epsilon] = ...
        dcHMCRunOnce(task.f,task.nTau,...
        task.x,task.g,task.E,task.r,task.epsilon);
    task.zTrials = task.zTrials + 1;
    
    % create new sample
    if mod(task.zTrials,task.thin) == 0
        task.zSamples = task.zSamples + 1;
        task.samples(task.zSamples,:) = task.x;
        task.L(task.zSamples) = -task.E;
    end
end
result = task;
end

