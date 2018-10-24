function varargout=xyzlabels(ax,xl,yl,zl,varargin)
% Creator: Joshua Yeh
% Date created: 2018-01-09
% 
%% DESCRIPTION
% The purpose of this script is to consolidate the xlabel and ylabel
% built-in MATLAB functions in order to create more concise coding.
% 
%% INPUT VARIABLES
% ax: axes handles in which the labels will be placed in (single col or row
% array)
% xl: xlabel string
% yl: ylabel string
% zl: zlabel string
% varargin:
%     "fontweight": define fontweight for both labels
%     "fontsize": define fontsize for both labels
% 
%% OUTPUT VARIABLES
% varargout{1}: handle to xlabel
% varargout{2}: handle to ylabel

% Check to see if xl and yl are character strings
if ~ischar(xl)||~ischar(yl)||~ischar(zl)
    error('xl, yl, and zl must be character strings');
end

% Parse input variables
narginchk(4,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('fontweight','normal',@(x) ischar(x));
params.addParameter('fontsize',ax.FontSize,@(x) isnumeric(x));
params.addParameter('fontname','Microsoft YaHei Light',@(x) ischar(x));
params.addParameter('interpreter','tex',@(x) strcmp(x,'tex')|strcmp(x,'latex'));
params.parse(varargin{:});

% extracted parsed parameters
fontweight=params.Results.fontweight;%fontweight of texts
fontsize=params.Results.fontsize;%fontsize of texts
fontname=params.Results.fontname;%fontname
interpreter=params.Results.interpreter;%text interpreter

for dum=1:length(ax)
    xh0=xlabel(ax(dum),xl,'fontsize',fontsize,'fontweight',fontweight,...
        'fontname',fontname,'interpreter',interpreter);
    yh0=ylabel(ax(dum),yl,'fontsize',fontsize,'fontweight',fontweight,...
        'fontname',fontname,'interpreter',interpreter);
    zh0=zlabel(ax(dum),yl,'fontsize',fontsize,'fontweight',fontweight,...
        'fontname',fontname,'interpreter',interpreter);
    set(xh0,'units','normalized');
    set(yh0,'units','normalized');
    set(zh0,'units','normalized');
end

% define output
varargout{1}=xh0;%xlabel handle
varargout{2}=yh0;%ylabel handle
varargout{3}=zh0;%zlabel handle
