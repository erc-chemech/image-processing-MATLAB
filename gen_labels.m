function [xticklabels,yticklabels]=gen_labels(x,y)
% Author: Joshua Yeh
% Date created: 2017/11/02
% 
%% DESCRIPTION
% This function generates user defined xlabels and ylabels, specified by x
% and y.
% 
%% INPUT VARIABLES
% x: tick labels for the x axis (double)
% y: tick labels for the y axis (double)s
%
%% OUTPUT VARIABLES
% xticklabels: x tick labels (cell of strings)
% yticklabels: y tick labels (cell of strings)
% 
%%
% xlabels
xticklabels={'1'};%dummy label
count=1;
for dum=x
    xticklabels(count)={num2str(dum)};
    count=count+1;
end
xticklabels=xticklabels';%make into a vertical array

%y labels
yticklabels={'1'};%dummy label
count=1;
for dum=y
    yticklabels(count)={num2str(dum)};
    count=count+1;
end
yticklabels=yticklabels';%make into a vertical array