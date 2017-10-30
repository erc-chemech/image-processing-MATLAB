function output=import_tiff_stack(file)
%% DESCRIPTION
% This function imports a stacked tiff file and stores in an ouput variable
% structure.
%
% INPUT VARIABLES
% file: names of the tiff file to be imported
disp(['Importing ',file]);
tiff_info=imfinfo(file);%tiff info
tiff_stack=imread(file,1);%read first image in the stack
tiff_stack=medfilt2(tiff_stack,[10 10]);% remove noise using medium filter
I_sum_z=sum(sum(tiff_stack));% integrate intensity for first image plane

%extract and apply median filter to each image in stack and concatenate
for dum=2:size(tiff_info,1)
    plane=medfilt2(imread(file,dum),[10 10]);
    tiff_stack=cat(3,tiff_stack,plane);
    I_sum_z=cat(1,I_sum_z,sum(sum(plane)));%integrate intensity for image plane
end
tiff_stack=double(tiff_stack);%convert to double (64-bit precision)
tiff_stack_sum=sum(tiff_stack,3);%sum all of the planes

%Store variables in the output variable structure
output.tiff_stack=tiff_stack;%filtered (via medfilt2) 3d image
output.I_sum_z=I_sum_z;%total intensity for each plane
output.tiff_stack_sum=tiff_stack_sum;%summed intensity through thickness
output.file=file;%filename of the imported tiff file
output.info=imfinfo(file);%Stores the metadata information
disp('Import finished');
