function displayResultsForTask(results,printDir)
if ~exist('printDir','var')
    printDir = [];
end

if isfield(results,'dataGE')
    printSimultaneousResults(results,printDir);
else
    printSingleResults(results,printDir);
end
end

function printSingleResults(results,printDir)
zSamples = results.zSamples;
L           = results.L(1:zSamples);
samples     = results.samples(1:zSamples,:);
sim         = results.sim;
beta        = results.beta;
roi         = results.roi;
data        = results.data;
sesname     = results.sesname;
noBVMaskFlag= results.noBVMaskFlag;

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

if ~isempty(printDir)
    fileName = fullfile(printDir,[sesname noBVString '_results.txt' ]);
    folderName = fileparts(fileName);
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end
    fid = fopen(fileName,'w+');
else
    fid = 1;
end
    
[~,medianPSample] = min(abs(L-median(L)));
[~,maxPSample] = max(L);

[noiseSamples,rhoSamples,deltaSamples,epsilonSamples,...
    omegaSamples,thetaSamples,fwhmSamples] = ...
    vectorToParameters(samples,sim);

f = [1/12 1/1.5];
R2 = zeros(1,zSamples);
R2filtered = zeros(1,zSamples);
dataFiltered = filterMap(data,roi,f,sim);

for z=1:zSamples
    [~,mriPattern] = ...
        computeODCModel(sim,squeeze(noiseSamples(z,:,:)),...
        rhoSamples(z),deltaSamples(z),...
        epsilonSamples(z),thetaSamples(z),...
        omegaSamples(z),beta,fwhmSamples(z));
    cc = corrcoef(data(roi==1),mriPattern(roi==1));
    R2(z) = cc(1,2).^2;
    
    mriPatternFiltered = filterMap(mriPattern,roi,f,sim);
    cc = corrcoef(dataFiltered(roi==1),mriPatternFiltered(roi==1));
    R2filtered(z) = cc(1,2).^2;
end

fwhm = mean(fwhmSamples);
fwhmstd = std(fwhmSamples);
fprintf(fid,'%s\n',sesname);
fprintf(fid,'fwhm = %.2f +- %.2f\n',fwhm,fwhmstd);
fprintf(fid,'beta = %.3f\n',results.beta);
fprintf(fid,'R^2 of median posterior sample = %.2f\n',R2(medianPSample));
fprintf(fid,'R^2 of maximum posterior sample = %.2f\n',R2(maxPSample));
fprintf(fid,'mean R^2 = %.2f\n',mean(R2));
fprintf(fid,'max R^2 = %.2f\n',max(R2));
fprintf(fid,'filtered R^2 of median posterior sample = %.2f\n',...
    R2filtered(medianPSample));
fprintf(fid,'filtered R^2 of maximum posterior sample = %.2f\n',...
    R2filtered(maxPSample));
fprintf(fid,'filtered mean R^2 = %.2f\n',...
    mean(R2filtered));
fprintf(fid,'filtered max R^2 = %.2f\n',...
    max(R2filtered));

if fid~=1
    fclose(fid);
end
end

function printSimultaneousResults(results,printDir)
zSamples = results.zSamples;
L           = results.L(1:zSamples);
samples     = results.samples(1:zSamples,:);
sim = results.sim;
betaGE = results.betaGE;
betaSE = results.betaSE;
roiGE = results.roiGE;
roiSE = results.roiSE;
dataGE        = results.dataGE;
dataSE        = results.dataSE;
noBVMaskFlag= results.noBVMaskFlag;

sesname = [results.sesnameGE '_' results.sesnameSE];

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

if ~isempty(printDir)
    fileName = fullfile(printDir,[sesname noBVString '_results.txt' ]);
    folderName = fileparts(fileName);
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end
    fid = fopen(fileName,'w+');
else
    fid = 1;
end

[~,medianPSample] = min(abs(L-median(L)));
[~,maxPSample] = max(L);

[noiseSamples,rhoSamples,deltaSamples,epsilonSamples,...
    omegaSamples,thetaSamples,fwhmSamplesGE,fwhmSamplesSE] = ...
    vectorToParametersSimultaneous(samples,sim);

R2GE = zeros(1,zSamples);
R2SE = zeros(1,zSamples);


f = [1/12 1/1.5];
R2GEfiltered = zeros(1,zSamples);
R2SEfiltered = zeros(1,zSamples);
dataGEFiltered = filterMap(dataGE,roiGE,f,sim);
dataSEFiltered = filterMap(dataSE,roiSE,f,sim);
for z=1:zSamples
    [~,mriPatternGE,mriPatternSE] = ...
        computeODCModelSimultaneous(sim,squeeze(noiseSamples(z,:,:)),...
        rhoSamples(z),deltaSamples(z),...
        epsilonSamples(z),thetaSamples(z),...
        omegaSamples(z),betaGE,betaSE,...
        fwhmSamplesGE(z),fwhmSamplesSE(z));
    ccGE = corrcoef(dataGE(roiGE==1),mriPatternGE(roiGE==1));
    ccSE = corrcoef(dataSE(roiSE==1),mriPatternSE(roiSE==1));
    R2GE(z) = ccGE(1,2).^2;
    R2SE(z) = ccSE(1,2).^2;
    
    mriPatternGEFiltered = filterMap(mriPatternGE,roiGE,f,sim);
    mriPatternSEFiltered = filterMap(mriPatternSE,roiSE,f,sim);

    ccGE = corrcoef(dataGEFiltered(roiGE==1)...
        ,mriPatternGEFiltered(roiGE==1));
    R2GEfiltered(z) = ccGE(1,2).^2;
    ccSE = corrcoef(dataSEFiltered(roiSE==1)...
        ,mriPatternSEFiltered(roiSE==1));
    R2SEfiltered(z) = ccSE(1,2).^2;

end

fwhmGE = mean(fwhmSamplesGE);
fwhmGEstd = std(fwhmSamplesGE);
fwhmSE = mean(fwhmSamplesSE);
fwhmSEstd = std(fwhmSamplesSE);

roi = (roiGE>0.5) & (roiSE>0.5);
c = corrcoef(dataGE(roi(:)),dataSE(roi(:)));

cf = corrcoef(dataGEFiltered(roi(:)),dataSEFiltered(roi(:)));

fprintf(fid,'%s\n',sesname);
fprintf(fid,'Correlation between GE and SE = %.2f\n',c(1,2));
fprintf(fid,'filtered correlation between GE and SE = %.2f\n',cf(1,2));

fprintf(fid,'GE (simultaneous)\n');
fprintf(fid,'fwhm = %.2f +- %.2f\n',fwhmGE,fwhmGEstd);
fprintf(fid,'R^2 of median posterior sample = %.2f\n',R2GE(medianPSample));
fprintf(fid,'R^2 of maximum posterior sample = %.2f\n',R2GE(maxPSample));
fprintf(fid,'mean R^2 = %.2f\n',mean(R2GE));
fprintf(fid,'max R^2 = %.2f\n',max(R2GE));
fprintf(fid,'filtered R^2 of median posterior sample = %.2f\n',...
    R2GEfiltered(medianPSample));
fprintf(fid,'filtered R^2 of maximum posterior sample = %.2f\n',...
    R2GEfiltered(maxPSample));
fprintf(fid,'filtered mean R^2 = %.2f\n',...
    mean(R2GEfiltered));
fprintf(fid,'filtered max R^2 = %.2f\n',...
    max(R2GEfiltered));

fprintf(fid,'SE (simultaneous)\n');
fprintf(fid,'fwhm = %.2f +- %.2f\n',fwhmSE,fwhmSEstd);
fprintf(fid,'R^2 of median posterior sample = %.2f\n',R2SE(medianPSample));
fprintf(fid,'R^2 of maximum posterior sample = %.2f\n',R2SE(maxPSample));
fprintf(fid,'mean R^2 = %.2f\n',mean(R2SE));
fprintf(fid,'max R^2 = %.2f\n',max(R2SE));
fprintf(fid,'filtered R^2 of median posterior sample = %.2f\n',...
    R2SEfiltered(medianPSample));
fprintf(fid,'filtered R^2 of maximum posterior sample = %.2f\n',...
    R2SEfiltered(maxPSample));
fprintf(fid,'filtered mean R^2 = %.2f\n',...
    mean(R2SEfiltered));
fprintf(fid,'filtered max R^2 = %.2f\n',...
    max(R2SEfiltered));

if fid~=1
    fclose(fid);
end
end

function imgf = filterMap(img,roi,f,sim)
res = [sim.vSize sim.vSize];

img(~roi) = mean(img(roi(:)==true));
kImg = fftshift(fft2(img));

% filter img
if ~isempty(f)
    N1 = size(kImg,1);
    N2 = size(kImg,2);
    L1 = N1 * res(1);
    L2 = N2 * res(2);
    dk1 = 1/L1;
    dk2 = 1/L2;
    k1Range = (-ceil((N1-1)/2):floor((N1-1)/2))*dk1;
    k2Range = (-ceil((N2-1)/2):floor((N2-1)/2))*dk2;
    [k1,k2] = ndgrid(k1Range,k2Range);
    kr = sqrt(k1.^2 + k2.^2);
    filterF = (kr>=f(1)) & (kr<=f(2));
    kImg = kImg.*filterF;
end
imgf = ifft2(ifftshift(kImg),'symmetric');
end