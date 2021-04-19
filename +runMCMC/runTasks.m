function tasks = runTasks(tasks,jobName)
maxTimeInHours = 1 * 12;
maxThin = 256;
multiFlag = false;
tasks = arrayfun(@(x) ...
    setfield(x,'jobName',jobName),tasks);
runMCMC.submitTasksForGuilliminCompiled(tasks,jobName,...
    maxTimeInHours,maxThin,multiFlag)
end
