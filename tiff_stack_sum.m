function T_obj=tiff_stack_sum(filein,fileout,varargin)
%% DESCRIPTION
% This function imports a tiff stack file and export a tiff summed image of
% the tiff file through the thickness or stack of the multi-page tiff. This
% function allows for the summed tiff file to be rotated (optional)
% 
%% SYNTAX
% tiff_stack_sum(filein,fileout)
% tiff_stack_sum(filein,fileout,r)
% tiff_stack_sum(filein,fileout,r,rotation)
% 
%% INPUT VARIABLES
% filein: the file that is imported (the multi-page tiff)
%
% r (varargin{1}): px width for median filter process
% rotation (varargin{2}): rotate the image counter-clockwise (deg.)
% 
%% OUTPUT VARIABLES
% T_obj: the Tiff object associated with the fileout tiff

switch nargin
    case 2
        r=10;
        rotation=0;
    case 3
        r=varargin{1};
        rotation=0;
    case 4
        r=varargin{1};
        rotation=varargin{2};
end

% Import tiff file
tiff1=import_tiff_stack(filein,r);

% Create Tiff object
T_obj=Tiff(fileout,'w');

% Define required tag fields
tags.ImageLength=tiff1.info(1).Height;
tags.ImageWidth=tiff1.info(1).Width;
tags.Photometric=1;
tags.BitsPerSample=tiff1.info(1).BitDepth;
tags.SamplesPerPixel=tiff1.info(1).SamplesPerPixel;
tags.RowsPerStrip=tiff1.info(1).RowsPerStrip;
tags.PlanarConfiguration=1;
tags.Software='Matlab';
T_obj.setTag(tags);

% rotate the image
output=imrotate(tiff1.tiff_stack_sum,rotation);

% resize image to defined height and width
output=output(1:tiff1.info(1).Height,1:tiff1.info(1).Width);

% Wrtie the tiff file
T_obj.write(uint16(output));
T_obj.close();%close connection to tiff file