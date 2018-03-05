function [out,A,B]=coord2image(X,Y,Z1,w,type)
% Author: Joshua Yeh
% Date created: 2018-21-2
%
% DESCRIPTION
% This script converts xyz coordinates into an image array.
% 
%% INPUT VARIABLES
% X: x coordinates
% 
% Y: y coordinates
% 
% Z1: z coordinates
% 
% w: the bin width or size of each pixel in the image array
% 
% type:
    % 'mean': take the mean of the datapoints within the same bin
    % 'none': uses the default behavior of summing all of the datapoints
% 
%% OUTPUT VARIABLES
% out: the image array
% 
% A: x meshgrid variable associated with the image array
% 
% B: y meshgrid variable associated with the image array
% 
%%

% This fcns converts coordinates of intensities into an image array
n1=min(X):w:max(X);
n2=min(Y):w:max(Y);
a1=linspace(min(X),max(X),numel(n1));
b1=linspace(min(Y),max(Y),numel(n2));

% Perform the binning procedure
ar=interp1(a1,1:numel(a1),X(:),'nearest');
br=interp1(b1,1:numel(b1),Y(:),'nearest');

% Convert coordinates into an image array
if strcmp(type,'mean')
    out=accumarray([br ar],Z1(:),[numel(n2) numel(n1)],@(x) mean(x));
elseif strcmp(type,'none')
    out=accumarray([br ar],Z1(:),[numel(n2) numel(n1)]);
elseif strcmp(type,'squared')
    out=accumarray([br ar],Z1(:),[numel(n2) numel(n1)],@(x) sum(x).^2);
end

% Create corresponding X and Y array associated withthe image array
[A,B]=meshgrid(a1,b1);
