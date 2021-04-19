function visualizeMCGoF(results, figDir)
samples     = results.samples;
data        = results.data;
roi         = results.roi;
beta        = results.beta;
L           = results.L;
sesname     = results.sesname;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;
increaseFactor = results.increaseFactor;
sim = results.sim;

nVoxels1 = sim.nVoxels1/increaseFactor;
nVoxels2 = sim.nVoxels2/increaseFactor;

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultFigureColor','White');
set(0,'DefaultAxesColor','None');

if exist('bwrLinCMap.mat','file')
    load('bwrLinCMap.mat');
else
    load('~/matlab/odcpsfMCMC/bwrLinCMap.mat');
end

[noiseSamples,rhoSamples,deltaSamples,epsilonSamples,...
    omegaSamples,thetaSamples,fwhmSamples] = ...
    vectorToParameters(samples(1:zSamples,:),sim);

[~,medianSample] = min(abs(L(1:zSamples)-median(L(1:zSamples))));
[~,bestSample] = max(L(1:zSamples));

R2 = zeros(1,zSamples);
for z=1:zSamples
    [~,mriPattern] = ...
        computeODCModel(sim,squeeze(noiseSamples(z,:,:)),...
        rhoSamples(z),deltaSamples(z),...
        epsilonSamples(z),thetaSamples(z),...
        omegaSamples(z),beta,fwhmSamples(z));
    cc = corrcoef(data(roi==1),mriPattern(roi==1));
    R2(z) = cc(1,2).^2;
end

noisePatternBest = squeeze(noiseSamples(bestSample,:,:));
[~,mriPatternBest] = ...
    computeODCModel(sim,noisePatternBest,...
    rhoSamples(bestSample),deltaSamples(bestSample),...
    epsilonSamples(bestSample),thetaSamples(bestSample),...
    omegaSamples(bestSample),beta,fwhmSamples(bestSample));

noisePatternMedian = squeeze(noiseSamples(medianSample,:,:));
[~,mriPatternMedian] = ...
    computeODCModel(sim,noisePatternMedian,...
    rhoSamples(medianSample),deltaSamples(medianSample),...
    epsilonSamples(medianSample),thetaSamples(medianSample),...
    omegaSamples(medianSample),beta,fwhmSamples(medianSample));

c = corrcoef(data(roi==1),mriPatternBest(roi==1));
RsqBest = c(1,2).^2;
c = corrcoef(data(roi==1),mriPatternMedian(roi==1));
RsqMedian = c(1,2).^2;
meanRsq = mean(R2);
maxRsq = max(R2);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

dataStr = regexprep([sesname noBVString],'_',' ');

f = figure;
subplot(3,6,[1 2 7 8]);
imagesc(data(1:nVoxels1,1:nVoxels2)'.*roi(1:nVoxels1,1:nVoxels2)',...
    [-beta beta]);
axis image;axis off;title(['data(' dataStr ')']);

subplot(3,6,[3 4 9 10]);
imagesc(mriPatternBest(1:nVoxels1,1:nVoxels2)'.*...
    roi(1:nVoxels1,1:nVoxels2)',[-beta beta]);
axis image;axis off;title('best sample');

subplot(3,6,[5 6 11 12]);
imagesc(mriPatternMedian(1:nVoxels1,1:nVoxels2)'.*...
    roi(1:nVoxels1,1:nVoxels2)',[-beta beta]);
axis image;axis off;title('median sample');

subplot(3,6,[15 16]);
line([-beta beta],[-beta beta],'Color','Black');
hold on;
plot(data(roi==1),mriPatternBest(roi==1),'.');
axis([-beta beta -beta beta]);
axis square; box on; grid on;
xlabel('measured responses');
ylabel('modeled responses');
title({['sample R^2 = ' sprintf('%.2f',RsqBest)],[' (mean/max=' ...
    sprintf('%.2f',meanRsq) '/' sprintf('%.2f',maxRsq) ')']});

subplot(3,6,[17 18]);
line([-beta beta],[-beta beta],'Color','Black');
hold on;
plot(data(roi==1),mriPatternMedian(roi==1),'.');
axis([-beta beta -beta beta]);
axis square; box on; grid on;
xlabel('measured responses');
ylabel('modeled responses');
title({['sample R^2=' sprintf('%.2f',RsqMedian)],['(mean/max=' ...
    sprintf('%.2f',meanRsq) '/' sprintf('%.2f',maxRsq) ')']});

colormap(bwrLinCMap);

dcprintfig([sesname noBVString '_MC3_GoF' ],24,18,figDir);
close(f);
end
