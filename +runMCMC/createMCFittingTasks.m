function tasks = createMCFittingTasks(subjectList,noBVMaskFlag,...
    removeMean,limitVector,nsamples,thin)
randSeed = 21;

seslistS1 = {'s1_ge',...
    's1_mostrepro_ge',...
    's1_se',...
    's1_mostrepro_se'};

seslistS2 = {'s2_ge',...
    's2_mostrepro_ge',...
    's2_se',...
    's2_mostrepro_se'};

seslist = {};
for zSubj=1:length(subjectList)
switch subjectList{zSubj}
    case 's1'
        seslist = cat(2,seslist,seslistS1);
    case 's2'
        seslist = cat(2,seslist,seslistS2);
    case 'mostrepro'
        seslist = ...
            seslist(cellfun(@(x) ~isempty(x),...
            strfind(seslist,'mostrepro')));
    case 'alldays'
        seslist = ...
            seslist(cellfun(@(x) ...
            strcmp(x(4:5),'ge')||strcmp(x(4:5),'se')...
            ,seslist));
    case 'allaverages'
        seslist = ...
            seslist(cellfun(@(x) isnan(str2double(x(4:5))), ...
            seslist));
    otherwise
        seslist = cat(2,seslist,subjectList{zSubj});
end
end

if ~exist('noBVMaskFlag','var')
    noBVMaskFlag = [0 1];
end
nNoBVMaskFlags = length(noBVMaskFlag); 
tasks = struct;
for zSes = 1:length(seslist)
    for zNoBVMaskFlag = 1:nNoBVMaskFlags
    zTask = (zSes-1)*nNoBVMaskFlags+zNoBVMaskFlag;
        [tasks(zTask).data, tasks(zTask).roi, tasks(zTask).beta, ...
            tasks(zTask).noiseLevel,tasks(zTask).filterCutoffs] = ...
        loadDataForFitting(seslist{zSes},noBVMaskFlag(zNoBVMaskFlag),removeMean);
        tasks(zTask).noBVMaskFlag = noBVMaskFlag(zNoBVMaskFlag);
        tasks(zTask).limitVector = limitVector;
        tasks(zTask).removeMean = removeMean;
        tasks(zTask).nsamples = nsamples;
        tasks(zTask).thin = thin;
        tasks(zTask).randSeed = randSeed;
        tasks(zTask).sesname = seslist{zSes};
        tasks(zTask).taskID = zTask;
    end
end
end
