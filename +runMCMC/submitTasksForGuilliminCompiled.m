function submitTasksForGuilliminCompiled(tasks,jobName,...
    maxTimeInHours,maxThin,multiFlag)

nTasks = length(tasks);

% create directory
mkdir(jobName);

% save data
for z = 1:length(tasks)
    task = tasks(z); %#ok<NASGU>
    save(fullfile(jobName,num2str(z)),'task','-v7.3');
end

% create submission script
if ~exist('multiFlag','var')
    multiFlag = false;
end
writeSubmissionScript(jobName,nTasks,maxTimeInHours,maxThin,multiFlag);

% transfer job folder (do not overwrite)
system(['rsync -az --ignore-existing ' ...
    jobName ' dchaimow@guillimin-p2.hpc.mcgill.ca:' ...
    fullfile('/sb/project/ngt-630-aa/denis/odcpsf/jobs/')]);

% but do overwrite submission script
system(['scp ' fullfile(jobName,'myBatchFile.sh') ...
    ' dchaimow@guillimin-p2.hpc.mcgill.ca:' ...
    fullfile('/sb/project/ngt-630-aa/denis/odcpsf/jobs/',...
    jobName,'myBatchFile.sh')]);

% remove local job folder
rmdir(jobName,'s');

% submit job
disp(['submitting job, estimated running time: ' ...
        datestr((maxTimeInHours*3600+3600)/(24*60*60),'DD:HH:MM:SS')]);
system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca "qsub ' ...
    fullfile('/sb/project/ngt-630-aa/denis/odcpsf/jobs/',jobName,...
    'myBatchFile.sh') '"']);
end

function writeSubmissionScript(...
    jobName,nTasks,maxTimeInHours,maxThin,multiFlag)

codeFolder = 'odcpsfMCMC';

wallTimeInSec = maxTimeInHours*3600+3600;
timeStr = datestr(wallTimeInSec/86400,'DD:HH:MM:SS');
if multiFlag
    runCommand = 'run_odcpsfmcmcMultiThreaded.sh';
else
    runCommand = 'run_odcpsfmcmc.sh';
end    

str0 =  '#!/bin/bash';
str1 =  '#PBS -A ngt-630-aa';
str2 = ['#PBS -N ' jobName];
if multiFlag
    str3 =  '#PBS -l nodes=1:ppn=8';
else
    str3 =  '#PBS -l nodes=1:ppn=1';
end
str4 =  '#PBS -l pmem=3gb';
str5 = ['#PBS -l walltime=' timeStr];
str6 = ['#PBS -t 1-' num2str(nTasks)];
str7 = ['#PBS -o /sb/project/ngt-630-aa/denis/odcpsf/jobs/' jobName '/' ...
    jobName '.out']; 
str8 = ['#PBS -e /sb/project/ngt-630-aa/denis/odcpsf/jobs/' jobName '/' ...
    jobName '.err']; 
str9 =  '';
str10=  'export MATLAB=/software/applications/matlab-2012a';
str11=  'export PATH=$PATH:$MATLAB/bin';
str12=  ['cd /sb/project/ngt-630-aa/denis/odcpsf/jobs/' jobName];
str13=  ['$HOME/matlab/' codeFolder '/' runCommand ' ' ...
    '$MATLAB $PBS_ARRAYID ' ...
    num2str(maxTimeInHours) ' ' num2str(maxThin) ...
    '  >> output_${PBS_ARRAYID}.out 2>&1'];
str = [str0 '\n' str1 '\n' str2 '\n' str3 '\n' str4 '\n' str5 '\n' str6 ...
    '\n' str7 '\n' str8 '\n' str9 '\n' str10 '\n' str11 '\n' str12 '\n' ...
    str13 '\n'];

fid = fopen(fullfile(jobName,'myBatchFile.sh'),'w');
fprintf(fid, str);
fclose(fid);
end
