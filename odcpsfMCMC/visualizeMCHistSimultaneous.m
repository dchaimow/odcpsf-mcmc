function visualizeMCHistSimultaneous(results, figDir)
binSize = 0.05;

samples     = results.samples;
noBVMaskFlag= results.noBVMaskFlag;
zSamples    = results.zSamples;

sesname = [results.sesnameGE '_' results.sesnameSE];

sim = results.sim;

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultFigureColor','White');
set(0,'DefaultAxesColor','None');

if exist('bwrLinCMap.mat','file')
    load('bwrLinCMap.mat');
else
    load('~/matlab/odcpsfMCMC/bwrLinCMap.mat');
end

[~,~,~,~,...
    ~,~,fwhmSamplesGE,fwhmSamplesSE] = ...
    vectorToParametersSimultaneous(samples(1:zSamples,:),sim);

if noBVMaskFlag
    noBVString = '';
else
    noBVString = '_bv_removed';
end

d = (fwhmSamplesGE-mean(fwhmSamplesGE)).^2+...
    (fwhmSamplesSE-mean(fwhmSamplesSE)).^2;
[~,meanIndx] = min(d);

f = figure;
subplot(2,2,2);
scatter(fwhmSamplesGE,fwhmSamplesSE,'MarkerEdgeColor',[0 0 0]);
hold on;
plot(fwhmSamplesGE(meanIndx),fwhmSamplesSE(meanIndx),'o',...
    'MarkerEdgeColor',[0 0.7 0]);
axis equal;
axis([0 2 0 2]);
line([0 2],[0 2],'Color',[0 0 0]);

subplot(2,2,4);
x = 0+binSize/2:binSize:2-binSize/2;
[n,x] = hist(fwhmSamplesGE,x);
bar(x,n,1,'FaceColor',[0.2 0.2 0.2],'EdgeColor','w');
axis([0 2 0 max(n)]);
line([mean(fwhmSamplesGE) mean(fwhmSamplesGE)],[0 max(n)],...
    'Color',[0 0.7 0]);
box off;
set(gca,'box','off','YTick',[],'YColor','w'); 
xlabel(['fwhm_{GE} (' sprintf('%.2f',mean(fwhmSamplesGE)) ...
    '\pm ' sprintf('%.2f',std(fwhmSamplesGE)) ')']);

subplot(2,2,1);
x = 0+binSize/2:binSize:2-binSize/2;
[n,x] = hist(fwhmSamplesSE,x);
barh(x,n,1,'FaceColor',[0.2 0.2 0.2],'EdgeColor','w');
axis([0 max(n) 0 2 ]);
line([0 max(n)],[mean(fwhmSamplesSE) mean(fwhmSamplesSE)],...
    'Color',[0 0.7 0]);
box off;
set(gca,'box','off','XTick',[],'XColor','w'); 
ylabel(['fwhm_{SE} (' sprintf('%.2f',mean(fwhmSamplesSE)) ...
    '\pm ' sprintf('%.2f',std(fwhmSamplesSE)) ')']);


dcprintfig([sesname noBVString '_MC2_Hist'],9,9,figDir);
close(f);
end


