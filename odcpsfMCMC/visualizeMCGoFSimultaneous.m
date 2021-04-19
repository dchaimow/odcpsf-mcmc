function visualizeMCGoFSimultaneous(results, figDir)
samples     = results.samples;
dataGE        = results.dataGE;
dataSE        = results.dataSE;
roiGE         = results.roiGE;
roiSE         = results.roiSE;
betaGE        = results.betaGE;
betaSE        = results.betaSE;
L           = results.L;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;
increaseFactor = results.increaseFactor;
sim = results.sim;

sesname = [results.sesnameGE '_' results.sesnameSE];

nVoxels1 = sim.nVoxels1/increaseFactor;
nVoxels2 = sim.nVoxels2/increaseFactor;

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultFigureColor','White');
set(0,'DefaultAxesColor','None');

if exist('bwrLinCMap.mat','file')
    load('bwrLinCMap.mat');
else
    load('~/matlab/odcpsfMCMC/bwrLinCMap.mat');
end

[noiseSamples,rhoSamples,deltaSamples,epsilonSamples,...
    omegaSamples,thetaSamples,fwhmSamplesGE,fwhmSamplesSE] = ...
    vectorToParametersSimultaneous(samples(1:zSamples,:),sim);

[~,bestSample] = max(L(1:zSamples));
[~,medianSample] = min(abs(L(1:zSamples)-median(L(1:zSamples))));

noisePatternBest = squeeze(noiseSamples(bestSample,:,:));
[~,mriPatternBestGE,mriPatternBestSE] = ...
    computeODCModelSimultaneous(sim,noisePatternBest,...
    rhoSamples(bestSample),deltaSamples(bestSample),...
    epsilonSamples(bestSample),thetaSamples(bestSample),...
    omegaSamples(bestSample),betaGE,betaSE,...
    fwhmSamplesGE(bestSample),fwhmSamplesSE(bestSample));

noisePatternMedian = squeeze(noiseSamples(medianSample,:,:));
[~,mriPatternMedianGE,mriPatternMedianSE] = ...
    computeODCModelSimultaneous(sim,noisePatternMedian,...
    rhoSamples(medianSample),deltaSamples(medianSample),...
    epsilonSamples(medianSample),thetaSamples(medianSample),...
    omegaSamples(medianSample),betaGE,betaSE,...
    fwhmSamplesGE(medianSample),fwhmSamplesSE(medianSample));

R2GE = zeros(1,zSamples);
R2SE = zeros(1,zSamples);
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
end

c = corrcoef(dataGE(roiGE==1),mriPatternBestGE(roiGE==1));
RsqBestGE = c(1,2).^2;
c = corrcoef(dataGE(roiGE==1),mriPatternMedianGE(roiGE==1));
RsqMedianGE = c(1,2).^2;

c = corrcoef(dataSE(roiSE==1),mriPatternBestSE(roiSE==1));
RsqBestSE = c(1,2).^2;
c = corrcoef(dataSE(roiSE==1),mriPatternMedianSE(roiSE==1));
RsqMedianSE = c(1,2).^2;

meanRsqGE = mean(R2GE);
meanRsqSE = mean(R2SE);
maxRsqGE = max(R2GE);
maxRsqSE = max(R2SE);



if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

f = figure;
subplot(6,6,[1 2 7 8]);
imagesc(dataGE(1:nVoxels1,1:nVoxels2)'.*roiGE(1:nVoxels1,1:nVoxels2)',...
    [-betaGE betaGE]);
axis image;axis off;title('data GE');

subplot(6,6,[3 4 9 10]);
imagesc(mriPatternBestGE(1:nVoxels1,1:nVoxels2)'.*...
    roiGE(1:nVoxels1,1:nVoxels2)',[-betaGE betaGE]);
axis image;axis off;title('best sample GE');

subplot(6,6,[5 6 11 12]);
imagesc(mriPatternMedianGE(1:nVoxels1,1:nVoxels2)'.*...
    roiGE(1:nVoxels1,1:nVoxels2)',[-betaGE betaGE]);
axis image;axis off;title('median sample GE');

subplot(6,6,[15 16]);
line([-betaGE betaGE],[-betaGE betaGE],'Color','Black');
hold on;
plot(dataGE(roiGE==1),mriPatternBestGE(roiGE==1),'.');
axis([-betaGE betaGE -betaGE betaGE]);
axis square; box on; grid on;
xlabel('measured GE');
ylabel('modeled GE');
title({['sample R^2=' sprintf('%.2f',RsqBestGE)],[' (mean/max=' ...
    sprintf('%.2f',meanRsqGE) '/' sprintf('%.2f',maxRsqGE) ')']});

subplot(6,6,[17 18]);
line([-betaGE betaGE],[-betaGE betaGE],'Color','Black');
hold on;
plot(dataGE(roiGE==1),mriPatternMedianGE(roiGE==1),'.');
axis([-betaGE betaGE -betaGE betaGE]);
axis square; box on; grid on;
xlabel('measured GE');
ylabel('modeled GE');
title({['sample R^2=' sprintf('%.2f',RsqMedianGE)],[' (mean/max=' ...
    sprintf('%.2f',meanRsqGE) '/' sprintf('%.2f',maxRsqGE) ')']});

subplot(6,6,[19 20 25 26]);
imagesc(dataSE(1:nVoxels1,1:nVoxels2)'.*roiSE(1:nVoxels1,1:nVoxels2)',...
    [-betaSE betaSE]);
axis image;axis off;title('data SE');

subplot(6,6,[21 22 27 28]);
imagesc(mriPatternBestSE(1:nVoxels1,1:nVoxels2)'.*...
    roiSE(1:nVoxels1,1:nVoxels2)',[-betaSE betaSE]);
axis image;axis off;title('best sample SE');

subplot(6,6,[23 24 29 30]);
imagesc(mriPatternMedianSE(1:nVoxels1,1:nVoxels2)'.*...
    roiSE(1:nVoxels1,1:nVoxels2)',[-betaSE betaSE]);
axis image;axis off;title('median sample SE');

subplot(6,6,[33 34]);
line([-betaSE betaSE],[-betaSE betaSE],'Color','Black');
hold on;
plot(dataSE(roiSE==1),mriPatternBestSE(roiSE==1),'.');
axis([-betaSE betaSE -betaSE betaSE]);
axis square; box on; grid on;
xlabel('measured SE');
ylabel('modeled SE');
title({['sample R^2=' sprintf('%.2f',RsqBestSE)],[' (mean/max=' ...
    sprintf('%.2f',meanRsqSE) '/' sprintf('%.2f',maxRsqSE) ')']});


subplot(6,6,[35 36]);
line([-betaSE betaSE],[-betaSE betaSE],'Color','Black');
hold on;
plot(dataSE(roiSE==1),mriPatternMedianSE(roiSE==1),'.');
axis([-betaSE betaSE -betaSE betaSE]);
axis square; box on; grid on;
xlabel('measured SE');
ylabel('modeled SE');
title({['sample R^2=' sprintf('%.2f',RsqMedianSE)],[' (mean/max=' ...
    sprintf('%.2f',meanRsqSE) '/' sprintf('%.2f',maxRsqSE) ')']});

colormap(bwrLinCMap);

dcprintfig([sesname noBVString '_MC3_GoF' ],32,48,figDir);
close(f);
end
