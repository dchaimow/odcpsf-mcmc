function visualizeMCSamples(results, figDir)
samples     = results.samples;
data        = results.data;
roi         = results.roi;
beta        = results.beta;
L           = results.L;
r           = results.r;
sesname     = results.sesname;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;
increaseFactor = results.increaseFactor;
sim = results.sim;

nVoxels1 = sim.nVoxels1/increaseFactor;
nVoxels2 = sim.nVoxels2/increaseFactor;
nSim1 = sim.N1/increaseFactor;
nSim2 = sim.N2/increaseFactor;

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

noisePattern = squeeze(noiseSamples(end,:,:));

[odcPattern,mriPattern] = computeODCModel(sim,noisePattern,...
    rhoSamples(end),deltaSamples(end),epsilonSamples(end),...
    thetaSamples(end),omegaSamples(end),beta,fwhmSamples(end));

roiUp = (dcupsample(roi,sim.upSampleFactor)>0.5);

f = figure;
subplot(4,6,[5 11]);
imagesc(noisePattern(1:nSim1,1:nSim2)'.*roiUp(1:nSim1,1:nSim2)',[-.1 .1]);
axis image;axis off;title('noise(end)');

subplot(4,6,[6 12]);
imagesc(odcPattern(1:nSim1,1:nSim2)'.*roiUp(1:nSim1,1:nSim2)',[-1 1]);
axis image;axis off;title('odc(end)');

subplot(4,6,[17 23]);
imagesc(mriPattern(1:nVoxels1,1:nVoxels2)'.*...
    roi(1:nVoxels1,1:nVoxels2)',[-beta beta]);
axis image;axis off;title('mri(end)');

subplot(4,6,[18 24]);
imagesc(data(1:nVoxels1,1:nVoxels2)'.*...
    roi(1:nVoxels1,1:nVoxels2)',[-beta beta]);
axis image;axis off;title('data');

subplot(4,6,[1 2]); plot(L(1:zSamples)); ...
    title(['logL (r=' num2str(r,2) ') thin=' num2str(results.thin)]);
    grid on;
subplot(4,6,[7 8]); plot(rhoSamples); ...
    title('rho'); grid on;
subplot(4,6,[13 14]); plot(deltaSamples); ...
    title('delta'); grid on;
subplot(4,6,[19 20]); plot(epsilonSamples); ...
    title('epsilon'); grid on;
subplot(4,6,[9 10]);  plot(omegaSamples); ...
    title('omega'); grid on;
subplot(4,6,[15 16]); plot(thetaSamples); ...
    title('theta'); grid on;
subplot(4,6,[21 22]); plot(fwhmSamples); ...
    title('fwhm'); grid on;

colormap(bwrLinCMap);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

dcprintfig([sesname noBVString '_MC1_Samples'],30,21,figDir);
close(f);
end
