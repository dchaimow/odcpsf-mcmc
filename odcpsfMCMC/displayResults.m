function displayResults(resultsFolder, writeFlag)
if exist('writeFlag','var') && writeFlag
    printDir = resultsFolder;
else
    printDir = [];
end
subjectList = {'s1','s2'};

for z=1:length(subjectList)
    load(fullfile(resultsFolder,[subjectList{z} '_ge']));
    displayResultsForTask(task,printDir);
    
    load(fullfile(resultsFolder,[subjectList{z} '_se']));
    displayResultsForTask(task,printDir);
    
    load(fullfile(resultsFolder,[subjectList{z} '_simultaneous']));
    displayResultsForTask(task,printDir);
end

end
