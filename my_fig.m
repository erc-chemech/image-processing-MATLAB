function f=my_fig(n,ax,varargin)
% Author: Joshua Yeh
% Date created: 2018/01/08
% 
%% DESCRIPTION
% This function creates a figure and outputs a figure structure that makes
% accessing the different axes easy. Note that this function is built ontop
% of the subtighplot function. Make sure that the subtightplot function is
% included in the path directory.
% 
%% INPUT VARIABLES
% n: figure number
% 
% ax: m element cell, where m defines the number of subplots. The column
% elements follows the same format as the elements defined in the MATLAB
% subplot function.
%   (Ex. ax{1}=[1,3,1] creates an
%   axes that is on a 1 by 3 grid and placed in the "1" position. For more
%   information, refer to the subplot function.
% 
% varargin: is parsed, where the options are as follows:
%   'gap': defines the gap between subaxes (used in subtightplot fcn)
%   'marg_h': margins in height in normalized units (used in subtightplot
%   fcn)
%   'marg_w': margins in width in normalized units (ised in subtightplot
%   fcn)
%   'fontsize': fontsize for all axes
% 
%% OUTPUT VARIABLES
% f: structure variables
%   f.s<axes number>: axes are numbered according to the order in which ax
%   is defined
%   f.f: handles to the figure
%
%% Parse input variables
narginchk(1,inf);

% If user doesn't specify the value for ax
if nargin==1
    ax={[1 1 1]};
end

params=inputParser;
params.CaseSensitive=false;
params.addParameter('gap',0.11,@(x) isnumeric(x));
params.addParameter('fontsize',18,@(x) isnumeric(x));
params.addParameter('fontname','Microsoft YaHei Light',@(x) ischar(x));
df=0;%flag for whether figure is a double axes plot
% if only 1 axes is inputed
if (numel(ax)==1&&isequal(ax{1},[1 1 1]))
    params.addParameter('marg_h',[0.17 0.1],@(x) length(x)==2&&isnumeric(x));
    params.addParameter('marg_w',[0.2 0.05],@(x) length(x)==2&&isnumeric(x));

% or double axes (2 axes stacked ontop of each other)
elseif (numel(ax)==2&&isequal(ax{1},[1 1 1]))&&isequal(ax{2},[1 1 1])
    params.addParameter('marg_h',[0.17 0.1],@(x) length(x)==2&&isnumeric(x));
    params.addParameter('marg_w',[0.17 0.17],@(x) length(x)==2&&isnumeric(x));
    df=1;
    
else%otherwise for multiple axes
    params.addParameter('marg_h',[0.17 0.1],@(x) length(x)==2&&isnumeric(x));
    params.addParameter('marg_w',[0.15 0.05],@(x) length(x)==2&&isnumeric(x));
end
params.parse(varargin{:});


%% Create figure object

% Extract values from input variables
gap=params.Results.gap;%gap between subplots
marg_h=params.Results.marg_h;%height margin used in subtightplot
marg_w=params.Results.marg_w;%width margin used in subtightplot
fontsize=params.Results.fontsize;%fontsize for all axes
fontname=params.Results.fontname;

% Define figure and store it in struct variable
f.f=figure(n); clf(f.f);

% Create subplots based on the ax input var
for dum=1:length(ax)
    fn=['s',num2str(dum)];%define fieldname for subplot
    
    if length(ax{dum})>4
        ax{dum}=ax{dum}(1:4);
    end
    
    % Create subplot
    f.(fn)=subtightplot(ax{dum}(1),ax{dum}(2),ax{dum}(3:end),gap,marg_h,marg_w);
    set(f.(fn),'nextplot','add','fontsize',fontsize,'fontname',fontname);
    
end

% If the plot is a double axes plot
if df==1
    f.s2=copyobj(f.s1,f.f);
    set(f.s1,'box','off');
    set(f.s2,'box','off','yaxislocation','right','xaxislocation','top',...
        'xticklabel',[],'color','none');
end
