function visualizeMCDependenciesSimultaneous(results,figDir)
samples     = results.samples;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;
sesname = [results.sesnameGE '_' results.sesnameSE];

sim = results.sim;

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultFigureColor','White');
set(0,'DefaultAxesColor','None');

if exist('bwrLinCMap.mat','file')
    load('bwrLinCMap.mat');
else
    load('~/matlab/odcpsfMCMC/bwrLinCMap.mat');
end

[~,rhoSamples,deltaSamples,epsilonSamples,...
    omegaSamples,~,fwhmSamplesGE,fwhmSamplesSE] = ...
    vectorToParametersSimultaneous(samples(1:zSamples,:),sim);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

f = figure;
subplot(4,2,1);
scatter(fwhmSamplesGE,omegaSamples);
xlabel('fwhm^{psf} GE [mm]');
ylabel('\omega');
grid on;
box on;
title('\omega(fwhm^{psf} GE)');

subplot(4,2,2);
scatter(fwhmSamplesGE,rhoSamples);
xlabel('fwhm^{psf} GE [mm]');
ylabel('\rho');
grid on;
box on;
title('\rho(fwhm^{psf} GE)');

subplot(4,2,3);
scatter(fwhmSamplesGE,deltaSamples);
xlabel('fwhm^{psf} [mm]');
ylabel('\delta');
grid on;
box on;
title('\delta(fwhm^{psf} GE)');

subplot(4,2,4);
scatter(fwhmSamplesGE,epsilonSamples);
xlabel('fwhm^{psf} GE [mm]');
ylabel('\epsilon');
grid on;
box on;
title('\epsilon(fwhm^{psf} GE)');



subplot(4,2,5);
scatter(fwhmSamplesSE,omegaSamples);
xlabel('fwhm^{psf} SE [mm]');
ylabel('\omega');
grid on;
box on;
title('\omega(fwhm^{psf} SE)');

subplot(4,2,6);
scatter(fwhmSamplesSE,rhoSamples);
xlabel('fwhm^{psf} SE [mm]');
ylabel('\rho');
grid on;
box on;
title('\rho(fwhm^{psf} SE)');

subplot(4,2,7);
scatter(fwhmSamplesSE,deltaSamples);
xlabel('fwhm^{psf} [mm]');
ylabel('\delta');
grid on;
box on;
title('\delta(fwhm^{psf} SE)');

subplot(4,2,8);
scatter(fwhmSamplesSE,epsilonSamples);
xlabel('fwhm^{psf} SE [mm]');
ylabel('\epsilon');
grid on;
box on;
title('\epsilon(fwhm^{psf} SE)');

dcprintfig([sesname noBVString '_MC4_Depend'],18,24,figDir);
close(f);
end


