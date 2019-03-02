function rgb_hist(I,a)
% Joshua Yeh
% Date: 190301
% 
%% DESCRIPTION
% The purpose of this fcn is to plot the rgb histogram of an rgb image
% array.
% 
%% INPUT VARIABLES
% I: image array to plot
% 
% a:axes
% 
%%

if nargin==1
    figure; axes;
elseif nargin==2
    axes(a);
end

hold on;
histogram(I(:,:,1),'facecolor','r');
histogram(I(:,:,2),'facecolor','g');
histogram(I(:,:,3),'facecolor','b');
xylabels(gca,'rgb intensity (a.u.)','counts');
center_axes(gca);
