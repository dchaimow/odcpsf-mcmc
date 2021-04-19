function ol = dcMADOutliers(data,crit)
% see:
% "Detecting outliers: Do not use standard deviation around the mean, 
% use absolute deviation around the median"
% Christophe Leys,Christophe Ley,Olivier Klein, Philippe Bernard,
% Laurent Licata
% Journal of experimental social psychology
% use 2.5 as default
maddata = mad(data,1)*1.4826;
ol = abs(data-nanmedian(data))>maddata*crit;
ol(isnan(data)) = 1;
ol(isinf(data)) = 1;
end
