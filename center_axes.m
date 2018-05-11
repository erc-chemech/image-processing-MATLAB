function center_axes(ax)
% Author: Joshua Yeh
% Date created: 2018/05/11
% 
%% DESCRIPTION
% This function accepts an axes handle and recenters the axes in the figure
% window.
% 
%% INPUT VARIABLES
% ax: axes handles
% 
%%
outerpos = ax.OuterPosition;%outer position of the axes
ti = ax.TightInset; %margins of text labels
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
