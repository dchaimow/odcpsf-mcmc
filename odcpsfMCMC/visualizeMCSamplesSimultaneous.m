function visualizeMCSamplesSimultaneous(results, figDir)
samples     = results.samples;
dataGE        = results.dataGE;
dataSE        = results.dataSE;
roiGE         = results.roiGE;
roiSE         = results.roiSE;
betaGE        = results.betaGE;
betaSE        = results.betaSE;
L           = results.L;
r           = results.r;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;
increaseFactor = results.increaseFactor;
sim = results.sim;

sesname = [results.sesnameGE '_' results.sesnameSE];


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
    omegaSamples,thetaSamples,fwhmSamplesGE,...
    fwhmSamplesSE,] = ...
    vectorToParametersSimultaneous(samples(1:zSamples,:),sim);

noisePattern = squeeze(noiseSamples(end,:,:));

[odcPattern,mriPatternGE,mriPatternSE] = ...
    computeODCModelSimultaneous(sim,noisePattern,...
    rhoSamples(end),deltaSamples(end),epsilonSamples(end),...
    thetaSamples(end),omegaSamples(end),betaGE,betaSE,...
    fwhmSamplesGE(end),fwhmSamplesSE(end));

roiGEUp = (dcupsample(roiGE,sim.upSampleFactor)>0.5);
roiSEUp = (dcupsample(roiSE,sim.upSampleFactor)>0.5);
roiUp = roiGEUp | roiSEUp;


f = figure;
subplot(4,8,[3 11]);
imagesc(noisePattern(1:nSim1,1:nSim2)'.*roiUp(1:nSim1,1:nSim2)',[-.1 .1]);
axis image;axis off;title('noise(end)');
subplot(4,8,[4 12]);
imagesc(odcPattern(1:nSim1,1:nSim2)'.*roiUp(1:nSim1,1:nSim2)',[-1 1]);
axis image;axis off;title('odc(end)');
subplot(4,8,[5 13]);
imagesc(mriPatternGE(1:nVoxels1,1:nVoxels2)'.*...
    roiGE(1:nVoxels1,1:nVoxels2)',[-betaGE betaGE]);
axis image;axis off;title('modelDataGE(end)');
subplot(4,8,[6 14]);
imagesc(dataGE(1:nVoxels1,1:nVoxels2)'.*...
    roiGE(1:nVoxels1,1:nVoxels2)',[-betaGE betaGE]);
axis image;axis off;title('dataGE');
subplot(4,8,[7 15]);
imagesc(mriPatternSE(1:nVoxels1,1:nVoxels2)'.*...
    roiSE(1:nVoxels1,1:nVoxels2)',[-betaSE betaSE]);
axis image;axis off;title('modelDataSE(end)');
subplot(4,8,[8 16]);
imagesc(dataSE(1:nVoxels1,1:nVoxels2)'.*...
    roiSE(1:nVoxels1,1:nVoxels2)',[-betaSE betaSE]);
axis image;axis off;title('dataSE');
subplot(4,8,[1 2]); plot(L(1:zSamples)); ...
    title(['logL (r=' num2str(r,2) ')']); grid on;
subplot(4,8,[9 10]); plot(rhoSamples); title('rho'); grid on;
subplot(4,8,[17 18]); plot(deltaSamples); title('delta'); grid on;
subplot(4,8,[19 20]); plot(epsilonSamples); title('epsilon'); grid on;
subplot(4,8,[25 26]);  plot(omegaSamples); title('omega'); grid on;
subplot(4,8,[27 28]); plot(thetaSamples); title('theta'); grid on;
subplot(4,8,[21 22]); plot(fwhmSamplesGE); title('fwhmGE'); grid on;
subplot(4,8,[23 24]); plot(fwhmSamplesSE); title('fwhmSE'); grid on;

colormap(bwrLinCMap);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

dcprintfig([sesname noBVString '_MC1_Samples'],24,18,figDir);
close(f);
end
