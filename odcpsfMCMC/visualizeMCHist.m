function visualizeMCHist(results, figDir)
binSize = 0.05;

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

[~,~,~,~,...
    ~,~,fwhmSamples] = ...
    vectorToParameters(samples(1:zSamples,:),sim);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end
dataSetString = regexprep([sesname noBVString],'_',' ');

f = figure;
x = 0+binSize/2:binSize:2-binSize/2;

[n,x] = hist(fwhmSamples,x);
bar(x,n,1,'FaceColor',[0.2 0.2 0.2],'EdgeColor','w');
axis([0 2 0 max(n)]);
line([mean(fwhmSamples) mean(fwhmSamples)],[0 max(n)],...
    'Color',[0 0.7 0]);
box off;
set(gca,'box','off','YTick',[],'YColor','w'); 
title([dataSetString ' fwhm_{psf}: ' sprintf('%.2f',mean(fwhmSamples)) ...
    '\pm ' sprintf('%.2f',std(fwhmSamples))]);
xlabel('mm');
ylabel('p[fwhm_{psf}|data]');

dcprintfig([sesname noBVString '_MC2_Hist'],9,6,figDir);
close(f);
end
