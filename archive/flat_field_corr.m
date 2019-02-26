function output=flat_field_corr(raw,B,white)
% Author: Joshua Yeh
% Date created: 2017/10/23
% 
%% DESCRIPTION
% Performs a flat field correction on a target image.
% 
%% INPUT VARIABLES
% raw: raw image
% B: "black" or "dark" background image
% white: "white" or "bright" reference image
% 
%% OUTPUT VARIABLES
% output: corrected image
% 
%%

image_corr=@(image,black,white) (image-black)./(white-black);
output=(raw-B)./(white-B);