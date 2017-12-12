function tiff_stack_sum(filein,fileout,r)
%% DESCRIPTION
% This function imports a tiff stack file and export a tiff summed image of
% the tiff file through the thickness or stack of the multi-page tiff.
% 
%% INPUT VARIABLES
% filein: the file that is imported (the multi-page tiff)
% 
%% OUTPUT VARIABLES
% fileout: the file that is exported (summed stacked image)
% 
switch nargin
    case 1
        r=10;
    case 2
        r=varargin{1};
end

% Import tiff file
tiff1=import_tiff_stack(filein,r);

% Create Tiff object
t=Tiff(fileout,'w');

% Define required tag fields
tags.ImageLength=tiff1.info.Height;
tags.ImageWidth=tiff1.info.Width;
tags.Photometric=1;
tags.BitsPerSample=tiff1.info.BitDepth;
tags.SamplesPerPixel=tiff1.info.SamplesPerPixel;
tags.RowsPerStrip=tiff1.info.RowsPerStrip;
tags.PlanarConfiguration=1;
tags.Software='Matlab';
t.setTag(tags);

% Wrtie the tiff file
t.write(uint16(tiff1.tiff_stack_sum));
t.close();%close connection to tiff file