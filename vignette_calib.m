function vignette_calib(varargin)
% Author: Joshua Yeh
% Date created: 2017/10/25
% 
%% VARARGIN
% vignette_calib(filenames)
%   filenames: structure variable that contains a list of filename strings
%   (including paths and extensions) 'filenames' must follow the following
%   format: {<tif images to be corrected>;...<black reference>;<white
%   reference>}
% 
% vignette_calib(filenames,image_out)
%   filenames: structure variable list of filename strings (including
%   paths)
%   image_out: the name of the output .mat file in which the imported tiff
%   images will be stored in or if image_out already exists, an attempt is
%   made to load the pre-imported image .mat file
% 
% vignette_calib(filenames,image_out,N)
% 
%% DESCRIPTION:
% This function script accepts a list of images according to a predefined
% format and performs a flat-field correction and background subtraction.  
% Images are scaled on grey scale (0 to 255) and saturated pixels (>255) 
% are replaced with a pixel value of 255. The number of saturated pixels 
% is displayed in the command window.
% Plots of the 'Nth' raw image, the white reference image, and the corrected
% images are shown along with a greyscale intensity of the corrected image.
%%

% Figure out what the variable inputs are
filenames=varargin{1};%import list of calibration files
switch nargin
    case 1        
        image_out='imported_images.mat';%import default desig. output filename
        N=1;%default index value, will plot the first image listed in filenames
    case 2
        [fp,f,ext]=fileparts(varargin{2});%identify fileparts of desig. output filename
        if isempty(f)==1
            image_out='imported_images.mat';%import default desig. output filename
        else
            image_out=[fp,'/',f,ext];%import desig. output filename
        end
        N=1;%default index value, will plot the first image listed in filenames
    case 3        
        [fp,f,ext]=fileparts(varargin{2});%identify fileparts of desig. output filename
        if isempty(fp)
            image_out=[f,ext];%import desig. output filename
        else
            image_out=[fp,'/',f,ext];%import desig. output filename
        end
        
        N=varargin{3};% index value, will plot the Nth image listed in filenames
end

try
    %% Import test image from mat file
    disp('Attempting to load any preimported image mat files');
    load(image_out,'calib');
    disp('Attempt successful');
catch
    disp('Error in importing mat file or mat file not found');
    disp('Attempting to import stacked tiff images...(this can take awhile)');
    %% Import the multi-stack for a standard calibraction image (fixed conc.)

    %Create a dummy 1x1 structure, this allows for concatenation in loop
    calib.tiff_stack=[];%filtered (via medfilt2) 3d image
    calib.I_sum_z=[];%total intensity for each plane
    calib.tiff_stack_sum=[];%summed intensity of planes through intensity
    calib.file=[];%filename of the imported tiff file
    calib.info=[];
    for dum=1:length(filenames)
        standard_file=filenames{dum};
        output=import_tiff_stack(standard_file);
        calib(dum)=output;
    end
    disp('Calibration image imported');
    disp('Saving workspace...');
    save(image_out);
    disp(['Workspace saved to ''',image_out,'''']);
end

%% Perform image correction
% Calibration is based on Eq. 2.11.14 from the book, "Current Protocols in
% Cytometry"

disp('Performing flat field correction');

image=calib(N).tiff_stack_sum;%image to be corrected (integrated through thickness)
black=calib(end-1).tiff_stack_sum;%black reference image (integrated through thickness)
white=calib(end).tiff_stack_sum;%white reference image (integrated through thickness)

% Perform image correction
image_corr=flat_field_corr(image,black,white);

disp('Correction completed');

%% plot raw image and corrected image
disp('Plotting...');

% Figure creation
f1.f=figure(1); clf(figure(1));
f1.f.Position=[15 520 1200 420];
f1.s1=subplot_tight(1,3,1,[0.1 0.05]);
f1.s2=subplot_tight(1,3,2,[0.1 0.05]);
f1.s3=subplot_tight(1,3,3,[0.1 0.05]);
set(findall(f1.f,'type','axes'),'nextplot','add');

% Show raw image
[~,f1.p1]=contour(f1.s1,image,25,'fill','on');
axis(f1.s1,'tight'); colorbar(f1.s1);
xlabel(f1.s1,'x');
ylabel(f1.s1,'y');
title(f1.s1,'raw image');
caxis(f1.s1,[0 max(white(:))]);

% Show raw white reference image
[~,f1.p2]=contour(f1.s2,white,25,'fill','on');
axis(f1.s2,'tight'); colorbar(f1.s2);
xlabel(f1.s2,'x');
ylabel(f1.s2,'y');
title(f1.s2,'standard');

% Show corrected image on greyscale (0 to 255)
[~,f1.p3]=contour(f1.s3,image_corr,linspace(1,1.5*nanmedian(image_corr(:)),25),...
    'fill','on');
axis(f1.s3,'tight'); colorbar(f1.s3);
xlabel(f1.s3,'x');
ylabel(f1.s3,'y');
title(f1.s3,'corrected image');
caxis(f1.s3,[0 255]);

%% Plot histograms of the corrected image

% Figure creation
f2.f=figure(2); clf(figure(2));
f2.f.Position=[360   500   560   315];
f2.f.Color='w';
f2.s1=axes;
set(findall(f2.f,'type','axes'),'nextplot','add','box','on');

% Histogram (greyscale) of corrected image
f2.p1=histogram(f2.s1,image_corr,0:255);
set(f2.p1,'edgecolor','none','facecolor','k');
set(f2.s1,'xlim',[0 255],'yscale','log','fontsize',16);
set(f2.s1,'ytick',...
    logspace(log10(f2.s1.YLim(1)),log10(f2.s1.YLim(2)),log10(f2.s1.YLim(2))+1));

% Show mean and median of distribution
f2.p2a=plot(f2.s1,ones(1,3).*nanmean(image_corr(:)),[f2.s1.YLim 1],...
    'color','r');
f2.p2b=plot(f2.s1,ones(1,3).*nanmedian(image_corr(:)),[f2.s1.YLim 1],...
    'color','b');

% Plot formatting
set(findall(f2.f,'type','line'),'linestyle','--','linewidth',2);
xlabel(f2.s1,'grey value intensity');
ylabel(f2.s1,'counts');
title(f2.s1,'intensity histogram');
L=legend([f2.p2a f2.p2b],['mean: ',num2str(nanmean(image_corr(:)))],...
    ['median: ',num2str(nanmedian(image_corr(:)))]);
set(findall(f2.f,'type','text'),'fontsize',16,'fontweight','bold');

disp('Plotting completed');