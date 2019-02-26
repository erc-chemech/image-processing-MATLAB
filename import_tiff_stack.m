function output=import_tiff_stack(file,r,varargin)
%% DESCRIPTION
% This function imports a stacked tiff file and stores in an ouput variable
% structure. This function only supports greyscale or binary images!
%
%% INPUT VARIABLES
% file: name of the tiff file to be imported (string variable)
% 
% r: the pixe side length of the square area in which a
% medfilt2 will be applied (setting r to 1 is the same as not applying a
% medfilt2 operation)
% 
% varargin 'fieldname',<value>
        % 'skip': skip basic calculations of images
        % 'silence': turn off some commands that display during importing
% 
%% OUTPUT VARIABLES
% output: structure variable containing the raw image file and integrated
%   intensity through the thickness of the stack
% 
%%
% parse user inputs
narginchk(2,inf);%check number of inputs is correct
params=inputParser;
params.addParameter('skip',0,@(x) isnumeric(x)&(x==1||x==0));
params.addParameter('silence',0,@(x) isnumeric(x)&(x==1||x==0));
params.parse(varargin{:});

skip=params.Results.skip;
silence=params.Results.silence;

% parse file variable
[~,name,ext]=fileparts(file);

disp(['Importing ',name,ext]);
tiff_info=imfinfo(file);%tiff info
disp(['Detected ',num2str(size(tiff_info,1)),' page(s) in tif file']);
tiff_stack=imread(file,1);%read first image in the stack
tiff_stack=medfilt2(double(tiff_stack),[r r]);% remove noise using medium filter

if skip==0
    I_sum_z=sum(tiff_stack(:));% integrate intensity for first image plane
    tiff_stack_std=std(tiff_stack(:));%determine std for each plane
    tiff_stack_mean=mean(tiff_stack(:));%determine mean for each plane
    tiff_stack_median=median(tiff_stack(:));%determine median for each plane
    tiff_stack_q1=quantile(tiff_stack(:),0.25);%determine 1st quartile
    tiff_stack_q3=quantile(tiff_stack(:),0.75);%determine 3rd quartile
end

%extract and apply median filter to each image in stack and concatenate
for dum=2:size(tiff_info,1)
    plane=double(medfilt2(imread(file,dum),[r r]));
    tiff_stack=cat(3,tiff_stack,plane);
    if skip==0
        I_sum_z=cat(1,I_sum_z,sum(plane(:)));%integrate intensity for image plane
        tiff_stack_std=cat(1,tiff_stack_std,std(plane(:)));%det. std for each plane
        tiff_stack_mean=cat(1,tiff_stack_mean,mean(plane(:)));%det. mean for each plane
        tiff_stack_median=cat(1,tiff_stack_median,median(tiff_stack(:)));
        tiff_stack_q1=cat(1,tiff_stack_q1,quantile(tiff_stack(:),0.25));
        tiff_stack_q3=cat(1,tiff_stack_q3,quantile(tiff_stack(:),0.75));
    end
    if mod(dum,ceil(size(tiff_info,1)*0.25))==0&&silence==0
        disp([num2str(dum),' of ',num2str(size(tiff_info,1)),' imported...']);
    end
end
% tiff_stack=double(tiff_stack);%convert to double
tiff_stack_sum=sum(tiff_stack,3);%determine sum of each plane

%Store variables in the output variable structure
output.tiff_stack=tiff_stack;%filtered (via medfilt2) 3d image
if skip==0
    output.I_sum_z=I_sum_z;%total intensity for each plane
    output.tiff_stack_sum=tiff_stack_sum;%summed intensity through thickness
    output.tiff_stack_mean=tiff_stack_mean;%mean of each plane
    output.tiff_stack_median=tiff_stack_median;%median of each plane
    output.tiff_stack_q1=tiff_stack_q1;%1st quartile of each plane
    output.tiff_stack_q3=tiff_stack_q3;%1st quartile of each plane
    output.tiff_stack_std=tiff_stack_std;%std of each plane
end
output.file=file;%filename of the imported tiff file
output.info=imfinfo(file);%Stores the metadata information
disp('Import finished');
