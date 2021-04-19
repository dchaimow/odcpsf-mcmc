function visualizeResults(results,figParentDir)
if ~exist('figParentDir','var')
    figParentDir = fullfile(...
        '~/Documents/Projects/currentProjects/odcpsf/Figures/',...
        [datestr(now,'yyyy_mm_dd') '_MCMCResults']);
end

for z=1:length(results)
    if isfield(results(z),'jobName')
        figDir = fullfile(figParentDir,results(z).jobName);
    else
        figDir = figParentDir;
    end
    if isfield(results(z),'dataGE')
        visualizeMCHistSimultaneous(results(z),figDir);
        visualizeMCDependenciesSimultaneous(results(z),figDir);
        visualizeMCGoFSimultaneous(results(z),figDir);
        visualizeMCSamplesSimultaneous(results(z),figDir);
    else
        visualizeMCHist(results(z),figDir);
        visualizeMCDependencies(results(z),figDir);
        visualizeMCGoF(results(z),figDir);
        visualizeMCSamples(results(z),figDir);
    end
end
end
