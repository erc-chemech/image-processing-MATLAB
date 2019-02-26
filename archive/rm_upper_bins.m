function gray_image_out=rm_upper_bins(gray_image,thresh,thresh_area)
% Author: Joshua Yeh
% Date created: 2017/10/31
%
%% DESCRIPTION
% This function accepts a gray_image on a 0-255 scale, discretizes the
% image into 255 bins, and removes the the lowest bin (setting those
% indices to nans). The purpose of this is to replace the background with
% nan values. This funcion also cleans then image from spurious nans and
% non-nans.
%
%% INPUT VARIABLES
% gray_image: a 2d gray image that will be processed
% thresh: user-defined threshold for background subtraction
% thresh_area: user-defined threshold for defining the area cutoff of spurious
% pixels
%% OUTPUT VARIABLES
% gray_image_out: a 2d (double) gray image with the background values as nan
% 
%% Discretize the gray image and replace background with nans
[m,n]=size(gray_image);%get size of image
[~,E]=discretize(double(gray_image),255);%discretize into 255 bins

%find first bin edge corresponding to thresh
kk=find(abs(E-thresh)==min(abs(E-thresh)));

%find the indices that were grouped into the 1st bin
ii=find(gray_image(:)>E(kk));
gray_image2=double(gray_image);
gray_image2(ii)=nan;%relace with nans

%% Now we try to get rid of spurious nans

% find all of the nans and store as a logical array
gray_image3=isnan(gray_image2);

% Get rid of small holes or small areas of nans. Threshold area is
% determined by 0.2% of the total pixels contained in the image.
thresh2=ceil(thresh_area*m*n);
gray_image4=bwareaopen(gray_image3,thresh2);
ii=find(gray_image4==1);%find all of the nan indices
gray_image5=double(gray_image);
gray_image5(ii)=nan;

%% Now we try to get rid of spurious non-nans likely associated with background

% find all of the non-nans and store as a logical array
gray_image6=isnan(gray_image5)==0;

% Get rid of small holes or small areas of non-nans. Threshold area is
% determined by 0.2% of the total pixels contained in the image.
gray_image7=bwareaopen(gray_image6,thresh2);
ii=find(gray_image7==0);%find all of the nan indices
gray_image_out=gray_image5;
gray_image_out(ii)=nan;