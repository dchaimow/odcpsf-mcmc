function tasks = createMCFittingTasksSimultaneous(subjectList,noBVMaskFlag,...
    removeMean,limitVector,nsamples,thin)
randSeed = 21;

seslistGE = {};
seslistSE = {};
for zSubj=1:length(subjectList)
    switch subjectList{zSubj}
        case 's2'
            seslistGE = cat(2,seslistGE,'s2_ge');
            seslistGE = cat(2,seslistGE,'s2_mostrepro_ge');
            seslistSE = cat(2,seslistSE,'s2_se');
            seslistSE = cat(2,seslistSE,'s2_mostrepro_se');
        case 'md'
            seslistGE = cat(2,seslistGE,'md_ge');
            seslistGE = cat(2,seslistGE,'md_mostrepro_ge');
            seslistSE = cat(2,seslistSE,'md_se');
            seslistSE = cat(2,seslistSE,'md_mostrepro_se');
        case 's1'
            seslistGE = cat(2,seslistGE,'s1_ge');
            seslistGE = cat(2,seslistGE,'s1_mostrepro_ge');
            seslistSE = cat(2,seslistSE,'s1_se');
            seslistSE = cat(2,seslistSE,'s1_mostrepro_se');
        case 'mostrepro'
            seslistGE = ...
                seslistGE(cellfun(@(x) ~isempty(x),...
                strfind(seslistGE,'mostrepro')));
            seslistSE = ...
                seslistSE(cellfun(@(x) ~isempty(x),...
                strfind(seslistSE,'mostrepro')));
        case 'alldays'
            seslistGE = ...
                seslistGE(cellfun(@(x) ...
                strcmp(x(4:5),'ge')...
                ,seslistGE));
            seslistSE = ...
                seslistSE(cellfun(@(x) ...
                strcmp(x(4:5),'se')...
                ,seslistSE));
        case 'allaverages'
            % dont remove anything
        otherwise
            if strfind(subjectList{zSubj},'_ge');
                seslistGE = cat(2,seslistGE,subjectList{zSubj});
            elseif strfind(subjectList{zSubj},'_se');
                seslistSE = cat(2,seslistSE,subjectList{zSubj});
            end
    end
end

nSes = length(seslistGE);
tasks = struct;

if length(noBVMaskFlag) == 1
    noBVMaskFlag = [noBVMaskFlag(1) noBVMaskFlag(1)];
end

for zTask = 1:nSes
    [tasks(zTask).dataGE, tasks(zTask).roiGE, tasks(zTask).betaGE, ...
        tasks(zTask).noiseLevelGE,tasks(zTask).filterCutoffsGE] = ...
        loadDataForFitting(...
        seslistGE{zTask},noBVMaskFlag(1),removeMean);
    [tasks(zTask).dataSE, tasks(zTask).roiSE, tasks(zTask).betaSE, ...
        tasks(zTask).noiseLevelSE,tasks(zTask).filterCutoffsSE] = ...
        loadDataForFitting(...
        seslistSE{zTask},noBVMaskFlag(2),removeMean);
    tasks(zTask).noBVMaskFlag = noBVMaskFlag;
    tasks(zTask).removeMean = removeMean;
    tasks(zTask).limitVector = limitVector;
    tasks(zTask).nsamples = nsamples;
    tasks(zTask).thin = thin;
    tasks(zTask).randSeed = randSeed;
    tasks(zTask).sesnameGE = seslistGE{zTask};
    tasks(zTask).sesnameSE = seslistSE{zTask};
    tasks(zTask).taskID = zTask;
end
end
