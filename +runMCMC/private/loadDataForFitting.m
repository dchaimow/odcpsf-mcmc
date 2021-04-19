function [diffMap, roi, beta, noiseLevel, filterCutoffs] = ...
    loadDataForFitting(sesname, noBVMaskFlag, removeMean)
addpath('../..');
dataDir = fullfile('data');
subjName = sesname(1:2);
sesDataDir = fullfile(dataDir,sesname);
subjDataDir = fullfile(dataDir,subjName);

if strcmp(sesname(3:end),'_ge')
    load(fullfile(subjDataDir,'results',[subjName '_glm_GE']));
    rcbvalues = gercbvalues;
    rcbevalues = gercbevalues;
    filterCutoffs = geFilterCutoffs;
elseif strcmp(sesname(3:end),'_se')
    load(fullfile(subjDataDir,'results',[subjName '_glm_SE']));
    rcbvalues = sercbvalues;
    rcbevalues = sercbevalues;
    filterCutoffs = seFilterCutoffs;
elseif strcmp(sesname(3:end),'_mostrepro_ge')
    load(fullfile(subjDataDir,'results',[subjName '_mostrepro_glm_GE']));
    rcbvalues = gercbvalues;
    rcbevalues = gercbevalues; 
    filterCutoffs = geFilterCutoffs;
elseif strcmp(sesname(3:end),'_mostrepro_se')
    load(fullfile(subjDataDir,'results',[subjName '_mostrepro_glm_SE']));
    rcbvalues = sercbvalues;
    rcbevalues = sercbevalues;
    filterCutoffs = seFilterCutoffs;
else
    load(fullfile(sesDataDir,'results',[sesname '_glm_allruns']));
    [roi,~] = cbiReadNifti(...
        fullfile(sesDataDir,'anat',[sesname '_roi_allreg.hdr']));
end
roi = roi>0.5;

if ~exist('filterCutoffs','var')
    filterCutoffs = [];
end

load(fullfile(subjDataDir,'results',[subjName '_bv_mask']));

diffvalues = rcbvalues(2,:);
if exist('noBVMaskFlag','var') && noBVMaskFlag
    bvMask = zeros(size(bvMask)); %#ok<NODEF>
end

outliersInNoBV = tools.dcMADOutliers(diffvalues(~bvMask),2.5);
outlierMask = zeros(size(bvMask));
outlierMask(~bvMask) = outliersInNoBV;
mask = bvMask | outlierMask;
% apply bv mask and outlier mask to data
diffvalues(mask) = 0;
% apply bv mask and outlier mask to roi
roiOriginal = roi; % save original roi
roi(roi(:)==1) = roi(roi(:)==1) & ~mask;

% remove mean from data
if exist('removeMean','var') && removeMean
    diffvalues(~mask) = diffvalues(~mask) - mean(diffvalues(~mask));
end

diffMap = tools.dccutroi(tools.dcunflat(diffvalues(~mask),roi),...
    roiOriginal);
roi = tools.dccutroi(roi,roiOriginal);


% %% remove outlier voxels in differential map
% outliersInNotBvMasked = dcMADOutliers(rcbvalues(2,~bvMask),2.5);
% outliers = zeros(size(bvMask));
% outliers(~bvMask) = outliersInNotBvMasked;
% mask = bvMask | outliers;
% 
% maskMap = dccutroi(dcunflat(mask',roi),roi);
% roi = dccutroi(roi,roi).*(~maskMap);

% beta calculation as twice the median of mean responses
beta = 2 * median(rcbvalues(1,~mask));
% noise level as RMS of standard error of differential response
noiseLevel = sqrt(mean(rcbevalues(2,~mask).^2));

% pad data if size is odd:
if mod(size(diffMap,1),2)
    diffMap = cat(1,zeros(1,size(diffMap,2)),diffMap);
    roi = cat(1,zeros(1,size(roi,2)),roi);
end

if mod(size(diffMap,2),2)
    diffMap = cat(2,zeros(size(diffMap,1),1),diffMap);
    roi = cat(2,zeros(size(roi,1),1),roi);
end
end
