function runSingleMCMCTask(fileName,maxTimeInHours,maxThin)
startTime = clock;
trialsBetweenSave = 2000;
figureDir = fullfile(fileparts(fileName),'../figures');

if ~exist('maxTimeInHours','var')
    maxTimeInHours = Inf;
end

if ~exist('maxThin','var')
    maxThin = Inf;
end

vars = load(fileName,'task');
task = vars.task;

if ~isfield(task,'samples')
    task = initializeTask(task);
    disp('initializing and starting new task');
else
    disp(['continuing aborted task at thin=' num2str(task.thin) ...
        ' zSamples=' num2str(task.zSamples)]);
end
 
rng(task.randState)     % set random number generator
while true
    while task.zSamples < 2*task.nsamples 
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
        
        % save results
        if mod(task.zTrials,trialsBetweenSave)==0
            task.randState = rng;
            save([fileName '_buf'],'task'); % first write to buffer file
            system(['mv ' fileName '_buf.mat ' fileName '.mat']); 
            visualizeResults(task,figureDir);
            elapsedTime = etime(clock,startTime)/3600;
            if elapsedTime>maxTimeInHours
                disp(['aborting (elapsed time ' ...
                    num2str(elapsedTime,2) 'h' ...
                    ' > maxTime=' num2str(maxTimeInHours) 'h)'])
                return;
            end
        end
    end
    if task.thin >= maxThin
        disp('aborting (maximum thin level completed)');
        task.randState = rng;
        save([fileName '_buf'],'task'); % first write to buffer file
        system(['mv ' fileName '_buf.mat ' fileName '.mat']);
        visualizeResults(task,figureDir); 
        return;
    end
    % prepare next level
    task.thin = task.thin * 2;
    task.samples(1:task.nsamples,:) = task.samples(2:2:2*task.nsamples,:);
    task.samples(task.nsamples+1:end,:) = 0;
    task.L(1:task.nsamples) = task.L(2:2:2*task.nsamples);
    task.L(task.nsamples+1:end) = 0;
    task.zTrials = 0;
    task.zSamples = task.nsamples;
    disp(['continuing with thin level ' num2str(task.thin)]);
end
end
