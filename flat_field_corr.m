% This fcn accepts an image and applies a flat field and background correction
% to the image. The corrected image is outputted.
% Joshua Yeh, 17/10/23
% raw: raw image
% B: "black" or "dark" background image
% white: "white" or "bright" reference image
function output=flat_field_corr(raw,B,white);
image_corr=@(image,black,white) (image-black)./(white-black);
output=(raw-B)./(white-B);