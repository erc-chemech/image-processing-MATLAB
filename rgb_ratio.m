function [R_ratio,G_ratio,B_ratio]=rgb_ratio(rgb_I)
% Author: Joshua Yeh
% Date created: 2017/10/31
% 
%% DESCRIPTION
% This function accepts an rgb image and computes the rgb ratios based on
% Grassman's Law from colorimetry. The image is deconstructed into
% normalized red, green, and blue images.
%
%% INPUT VARIABLES
% rgb_I: 2d rgb image that will be deconstructed into the normalized red,
% green, and blue channels (Grassman's Law). It is assumed that:
% rgb_I(:,:,1) is the red channel
% rgb_I(:,:,2) is the green channel
% rgb_I(:,:,3) is the blue channel
% 
%% OUTPUT VARIABLES
% R_ratio: normalized red channel (R_ratio=r/(r+b+g))
% G_ratio: normalized green channel (G_ratio=g/(r+b+g))
% B_ratio: normalized blue channel (B_ratio=b/(r+b+g))
% 
%%
% convert input image into double array
rgb_I=double(rgb_I);

% Find the sum of all three chaneels of the image
tot=sum(rgb_I,3);
r=rgb_I(:,:,1);%red channel
g=rgb_I(:,:,2);%green channel
b=rgb_I(:,:,3);%blue channel
R_ratio=r./tot;%normalized red channel
G_ratio=g./tot;%normalized green channel
B_ratio=b./tot;%normalized blue channel