function output=import_tiff_stack(file)
%% DESCRIPTION
% This function imports a stacked tiff file and stores in an ouput variable
% structure.
%
%% INPUT VARIABLES
% file: names of the tiff file to be imported
% 
%% OUTPUT VARIABLES
% output: structure variable contain the raw image file and integrated
%   intensity through the thickness of the stack
% 
%%

disp(['Importing ',file]);
tiff_info=imfinfo(file);%tiff info
disp(['Detected ',num2str(size(tiff_info,1)),' page(s) in tiff file']);
tiff_stack=imread(file,1);%read first image in the stack
tiff_stack=double(medfilt2(tiff_stack,[10 10]));% remove noise using medium filter
I_sum_z=sum(tiff_stack(:));% integrate intensity for first image plane
tiff_stack_std=std(tiff_stack(:));%determine std for each plane
tiff_stack_mean=mean(tiff_stack(:));%determine mean for each plane

%extract and apply median filter to each image in stack and concatenate
for dum=2:size(tiff_info,1)    
    plane=double(medfilt2(imread(file,dum),[10 10]));
    tiff_stack=cat(3,tiff_stack,plane);
    I_sum_z=cat(1,I_sum_z,sum(plane(:)));%integrate intensity for image plane
    tiff_stack_std=cat(1,tiff_stack_std,std(plane(:)));%det. std for each plane
    tiff_stack_mean=cat(1,tiff_stack_mean,mean(plane(:)));%det. mean for each plane
    if mod(dum,ceil(size(tiff_info,1)*0.75))==0
        disp('75% imported...');
    elseif mod(dum,ceil(size(tiff_info,1)*0.5))==0
        disp('50% imported...');
    end
end
% tiff_stack=double(tiff_stack);%convert to double
tiff_stack_sum=sum(tiff_stack,3);%determine sum of each plane

%Store variables in the output variable structure
output.tiff_stack=tiff_stack;%filtered (via medfilt2) 3d image
output.I_sum_z=I_sum_z;%total intensity for each plane
output.tiff_stack_sum=tiff_stack_sum;%summed intensity through thickness
output.tiff_stack_mean=tiff_stack_mean;%mean of each plane
output.tiff_stack_std=tiff_stack_std;%std of each plane
output.file=file;%filename of the imported tiff file
output.info=imfinfo(file);%Stores the metadata information
disp('Import finished');
