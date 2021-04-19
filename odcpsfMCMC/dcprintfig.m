function dcprintfig(description,w,v,figuredir,dateFlag)
%DCPRINTFIG  ...
%   Requirements
%
%   Design
%
%   Interfaces
%
%   Discussion
%

% Copyright 2009 Denis Chaimow.
% Created by Denis Chaimow on 12-May-2009 09:45:23
% $Id$'

%TODO: Write documentation
%TODO: Implement first version
%figuredir = '~/Desktop/figures';
if ~exist('figuredir','var') || isempty(figuredir)
    [ST,~] = dbstack;
    figuredir = fullfile('~/figures',ST(2).name);
end

if exist('dateFlag','var') && dateFlag
    [p,n,~] = fileparts(figuredir);
    figuredir = fullfile(p,...
        [datestr(now,'yyyy_mm_dd') '_' n]);
end

if ~exist(figuredir,'dir')
    mkdir(figuredir);
end

if ~exist('w','var');
    w = 9;
end
if ~exist('v','var');
    v = 9;
end

orient tall;
set(gcf, 'PaperPosition', [0 0 w v]);
%print('-depsc2',[fullfile(figuredir,description)]);
%print('-depsc2',[fullfile(figuredir,description) '_tmp']);
%system(['/usr/local/bin/eps2eps ' ...
%    fullfile(figuredir,description) '_tmp.eps ' ...
%    fullfile(figuredir,description) '.eps']);
%delete([fullfile(figuredir,description) '_tmp.eps']);
print('-djpeg','-r300',fullfile(figuredir,description));
%print('-dtiff','-r600',fullfile(figuredir,description));
%print('-depsc2','-painters',fullfile(figuredir,description));
%saveas(gcf,fullfile(figuredir,description));

%print('-depsc2',[fullfile(figuredir,description) '_tmp']);


%plot2svg(fullfile(figuredir,description), gcf, 'png'); 
