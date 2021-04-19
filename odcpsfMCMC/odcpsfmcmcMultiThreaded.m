function odcpsfmcmcMultiThreaded(fileName, maxTimeInHours, maxThin)
% compile using: 
% mcc -mv -R -nodisplay -N odcpsfmcmcMultiThreaded.m
if ~exist('maxTimeInHours','var')
    maxTimeInHours = Inf;
else
    maxTimeInHours = str2double(maxTimeInHours);
end

if ~exist('maxThin','var')
    maxThin = Inf;
else
    maxThin = str2double(maxThin);
end
runSingleMCMCTask(fileName,maxTimeInHours,maxThin);

