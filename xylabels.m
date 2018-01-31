function varargout=xylabels(ax,xl,yl,varargin)
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
% varargin:
%     "fontweight": define fontweight for both labels
%     "fontsize": define fontsize for both labels
% 
%% OUTPUT VARIABLES
% varargout{1}: handle to xlabel
% varargout{2}: handle to ylabel

% Check to see if xl and yl are character strings
if ~ischar(xl)&&~ischar(yl)
    error('xl and yl must be character string');
end

% Parse input variables
narginchk(3,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('fontweight','bold',@(x) ischar(x));
params.addParameter('fontsize',ax.FontSize,@(x) isnumeric(x));
params.parse(varargin{:});

% extracted parsed parameters
fontweight=params.Results.fontweight;%fontweight of texts
fontsize=params.Results.fontsize;%fontsize of texts

for dum=1:length(ax)
    xh0=xlabel(ax(dum),xl,'fontsize',fontsize,'fontweight',fontweight);
    yh0=ylabel(ax(dum),yl,'fontsize',fontsize,'fontweight',fontweight);
    set(xh0,'units','normalized');
    set(yh0,'units','normalized');
    
    % Need to work around bug in MATLAB saving, where the x and y labels
    % properties are not preserved
%     xh=copyobj(xh0,ax(dum));
%     yh=copyobj(yh0,ax(dum));
%     linkprop([xh ax(dum).XLabel],'position','fontsize');
%     linkprop([yh ax(dum).YLabel],'position','fontsize');
    
%     xlabel(ax(dum),[]);
%     ylabel(ax(dum),[]);
end

% define output
varargout{1}=xh0;%xlabel handle
varargout{2}=yh0;%ylabel handle