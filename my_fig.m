function f=my_fig(n,ax,varargin)
% Author: Joshua Yeh
% Date created: 2017/01/08
% 
%% DESCRIPTION
% This function creates a figure and outputs a figure structure that makes
% accessing the different axes easy. Note that this function is built ontop
% of the subtighplot function. Make sure that the subtightplkot function is
% included in the path directory.
% 
%% INPUT VARIABLES
% n: figure number
% 
% ax: m x 1 cell, where m defines the number of subplots. The column
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
% 
%% OUTPUT VARIABLES
% f: structure variables
%   f.s<axes number>: axes are numbered according to the order in which ax
%   is defined
%   f.f: handles to the figure

% Parse input variables
narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('gap',0.11,@(x) isnumeric(x));
params.addParameter('marg_h',[0.1 0.1],@(x) length(x)==2&&sum(isnumeric(x))==2);
params.addParameter('marg_w',[0.1 0.1],@(x) length(x)==2&&sum(isnumeric(x))==2);

% Extract values from input variables
gap=params.Results.gap;%gap between subplots
marg_h=params.Results.marg_h;%height margin used in subtightplot
marg_w=params.Results.marg_w;%width margin used in subtightplot

% Define figure and store it in struct variable
f.f=figure(n);

% Create subplots based on the ax input var
for dum=size(ax,1)
    fn=['s',num2str(dum)];%define fieldname
    
    % Create subplot
    f.(fn)=subtightplot(ax{dum}(1),ax{dum}(2),ax{dum}(3:end),gap,marg_h,marg_w);
end
