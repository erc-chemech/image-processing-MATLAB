% Author: Joshua Yeh
% Date created: 2017/10/25
% Example script on how to perform flat-field correction

% Preamble
clear all; close all; clc;
addpath('example images');%where the example tif files are stored

% Example of how the filenames structure should look like when using the
% vignette_calib MATLAB function
images={'20171019_acsn25b_ext0.tif';...
    '20171019_acsn50b_ext0.tif';...
    '20171019_acsn75b_ext0.tif';...
    '20171019_acsn0b_ext0.tif';...
    '20171019_acsn100b_ext0.tif'};

% This is where the imported tiff images will be saved as. Also,
% vignette_calib will search if this imported_tiff exists before running
% the flat_field_corr function in order to reduce run time.
imported_tiff='imported_calibration_images.mat';

% Index to be corrected and plotted.
N=1;

% In this case, 'AC-SN-0-B/20171019_acsn0b_ext0.tif' is the black reference
% image and 'AC-SN-100-B/20171019_acsn100b_ext0.tif' is the white
% reference. All other tif files are the ones that will be corrected
% according to the index,"N". Since N=1, the first image listed in
% standard_files will be corrected and plotted.

% Now run the vignette_calib function
vignette_calib(images,imported_tiff,1);