function odcpsfmcmc(fileName, maxTimeInHours, maxThin)
% compile using: 
% mcc -mv -R -singleCompThread -R -nodisplay -N odcpsfmcmc.m
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


