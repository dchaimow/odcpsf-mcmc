function visualizeMCDependencies(results,figDir)
samples     = results.samples;
sesname     = results.sesname;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;

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
    omegaSamples,~,fwhmSamples] = ...
    vectorToParameters(samples(1:zSamples,:),sim);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

f = figure;
subplot(2,2,1);
scatter(fwhmSamples,omegaSamples);
xlabel('fwhm^{psf} [mm]');
ylabel('\omega');
grid on;
box on;
title('\omega(fwhm^{psf})');

subplot(2,2,2);
scatter(fwhmSamples,rhoSamples);
xlabel('fwhm^{psf} [mm]');
ylabel('\rho');
grid on;
box on;
title('\rho(fwhm^{psf})');

subplot(2,2,3);
scatter(fwhmSamples,deltaSamples);
xlabel('fwhm^{psf} [mm]');
ylabel('\delta');
grid on;
box on;
title('\delta(fwhm^{psf})');

subplot(2,2,4);
scatter(fwhmSamples,epsilonSamples);
xlabel('fwhm^{psf} [mm]');
ylabel('\epsilon');
grid on;
box on;
title('\epsilon(fwhm^{psf})');

dcprintfig([sesname noBVString '_MC4_Depend'],18,12,figDir);
close(f);
end


