function T_obj=create_tiff(fileout,I)
%% DESCRIPTION
% This function that accepts an image, I, and exports the image as a tiff
% file. This function is different than the saveas built-in MATLAB
% function, since the saveas function grabs an rgb frame from an image.
% This function exports a matrix array as a tiff file.
% 
%% SYNTAX
% create_tiff(fileout,I)
% 
%% INPUT VARIABLES
% fileout: the file name in which the tiff image will be exported into
% 
% I: the image that will be converted into a tiff file
% 
%% OUTPUT VARIABLES
% T_obj: the Tiff object associated with the fileout tiff

% Create Tiff object
T_obj=Tiff(fileout,'w');

% Define required tag fields
tags.ImageLength=size(I,1);
tags.ImageWidth=size(I,2);
tags.Photometric=1;
tags.BitsPerSample=16;
tags.SamplesPerPixel=1;
tags.RowsPerStrip=8;
tags.PlanarConfiguration=1;
tags.Software='Matlab';
T_obj.setTag(tags);



% Wrtie the tiff file
T_obj.write(I);
T_obj.close();%close connection to tiff file

disp('Image exported to tiff file!');