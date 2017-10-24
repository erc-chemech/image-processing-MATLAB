% Joshua Yeh
% Calibration vignetting

%%Preamble
clear all; clc; close all;
addpath('../');

try
    %% Import test image from mat file
    disp('Attempting to load any prexisiting prepocessed image mat files');
    load('calibration_images.mat');
    B=calib(1);%"black" or "dark" reference image
    test=calib(2);%Image to be corrected
    disp('Attempt successful');
catch
    disp('Error in importing mat file or mat file not found');
    disp('Attempting to import stacked tiff images...');
    %% Import the multi-stack for a standard calibraction image (fixed conc.)
    standard_files={'../AC-SN-0-B/20171019_acsn0b_ext0.tif';...
        '../AC-SN-25-B/20171019_acsn25b_ext0.tif';...
        '../AC-SN-50-B/20171019_acsn50b_ext0.tif';...
        '../AC-SN-75-B/20171019_acsn75b_ext0.tif';...
        '../AC-SN-100-B/20171019_acsn100b_ext0.tif'};
    %Create a dummy 1x1 structure, this allows for concatenation in loop
    calib.tiff_stack=[];%filtered (via medfilt2) 3d image
    calib.I_sum_z=[];%total intensity for each plane
    calib.tiff_stack_sum=[];%summed intensity of planes through intensity
    calib.file=[];%filename of the imported tiff file
    calib.info=[];
    for dum=1:length(standard_files)
        standard_file=standard_files{dum};
        output=import_tiff_stack(standard_file);
        calib(dum)=output;
    end
    disp('Calibration image imported');
    disp('Saving workspace...');
    save('calibration_images.mat');
    disp('Workspace saved to ''calibration_images.mat''');
end

%% Perform image correction
% Calibration is based on Eq. 2.11.14 from the book, "Current Protocols in
% Cytometry"

disp('Performing flat field correction');

% Define function for image correction
image_corr=@(image,black,white) (image-black)./(white-black);

image=test.tiff_stack_sum;
white=calib(end).tiff_stack_sum;

%for plotting purposes only
white_corr=white./nanmax(white).*255;
ii=find(white_corr>255);
white_corr(ii)=nan;

test_corr=image_corr(image,...
    B.tiff_stack_sum,white).*255;
ii=find(test_corr>255);%get rid of greyscale values greater than max (255)
test_corr(ii)=nan;

disp('Correction completed');

%% plot raw image and corrected image
disp('Plotting...');

f1.f=figure(1); clf(figure(1));
f1.f.Position=[15 520 1200 420];
f1.s1=subplot_tight(1,3,1,[0.1 0.05]);
f1.s2=subplot_tight(1,3,2,[0.1 0.05]);
f1.s3=subplot_tight(1,3,3,[0.1 0.05]);
set(findall(f1.f,'type','axes'),'nextplot','add');

[~,f1.p1]=contour(f1.s1,image,25,'fill','on');
axis(f1.s1,'tight'); colorbar(f1.s1);
xlabel(f1.s1,'x');
ylabel(f1.s1,'y');
title(f1.s1,'raw image');
caxis(f1.s1,[0 max(white(:))]);

[~,f1.p2]=contour(f1.s2,white,25,'fill','on');
axis(f1.s2,'tight'); colorbar(f1.s2);
xlabel(f1.s2,'x');
ylabel(f1.s2,'y');
title(f1.s2,'standard');

[~,f1.p3]=contour(f1.s3,test_corr,linspace(1,1.5*nanmedian(test_corr(:)),25),...
    'fill','on');
axis(f1.s3,'tight'); colorbar(f1.s3);
xlabel(f1.s3,'x');
ylabel(f1.s3,'y');
title(f1.s3,'corrected image');
caxis(f1.s3,[0 255]);

%% Plot histograms of the corrected image
f2.f=figure(2); clf(figure(2));
f2.f.Position=[360   500   560   315];
f2.s1=axes;
set(findall(f2.f,'type','axes'),'nextplot','add','box','on');

f2.p1=histogram(f2.s1,test_corr,0:255);
set(f2.p1,'edgecolor','none','facecolor','k');
set(f2.s1,'xlim',[0 255],'yscale','log','fontsize',16);
set(f2.s1,'ytick',...
    logspace(log10(f2.s1.YLim(1)),log10(f2.s1.YLim(2)),log10(f2.s1.YLim(2))+1));

f2.p2a=plot(f2.s1,ones(1,3).*nanmean(test_corr(:)),[f2.s1.YLim 1],...
    'color','r');
f2.p2b=plot(f2.s1,ones(1,3).*nanmedian(test_corr(:)),[f2.s1.YLim 1],...
    'color','b');
set(findall(f2.f,'type','line'),'linestyle','--','linewidth',2);
xlabel(f2.s1,'grey value intensity');
ylabel(f2.s1,'counts');
title(f2.s1,'intensity histogram');
L=legend([f2.p2a f2.p2b],['mean: ',num2str(nanmean(test_corr(:)))],...
    ['median: ',num2str(nanmedian(test_corr(:)))]);
set(findall(f2.f,'type','text'),'fontsize',16,'fontweight','bold');

disp('Plotting completed');