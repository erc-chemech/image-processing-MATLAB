function [mov,f1]=script_frame_corr(varargin)
% Author: Joshua Yeh
% Date created: 2017/11/03
% 
%% DESCRIPTION
% This script automates frame correction for frames extracted from a
% mechanophore video experiment. A figure(1) window is created showing the
% original image of the region of interest (ROI) and a contour plot of the
% detected activated mechanophores in the sample. Histograms prior to
% correction and after correction is shown, along with the threshold value
% used to separate the signal (activated mechanophore) from the background.
% The corected image is based on the normalized RGB ratio (Grassman's Law)
% for the blue channel of the input RGB video.
% 
%% INPUT VARIABLES
% Must contain at least 4 variables
% filename: name of the file including the path
% frames: a uint8 single col or row array containing indices of the frames
% to be analyzed
% ROI: specifies the region of interest (ROI) of the original frame that
% will be analyzed.
% ref_ROI: specifies the reference region of interest (ROI) of the original
% frame that will be used as the white reference standard. This is needed
% to perform color correction of the frame.
% thresh: the threshold used to separate the signal (activated
% mechanophore) from the background of the frame. This value should be
% around 0.34.
% flag: a flag indicating whether or not to record the automated process of
% the frame correction (enable: 1, disable: 0, default: 0)
% 
%% OUTPUT VARIABLES
% mov: structure variables containing the frames that were extracted fromt
% he video file for analysis
% f1: figure handle of the plots created by this function script
% 
%% Determine what variables were inputted
switch nargin
    case 4
        filename=varargin{1};
        frames=varargin{2};
        ROI=varargin{3};
        ref_ROI=varargin{4};
        thresh=0.34;
        flag=0;
    case 5
        filename=varargin{1};
        frames=varargin{2};
        ROI=varargin{3};
        ref_ROI=varargin{4};
        thresh=varargin{5};
        flag=0;
    case 6
        filename=varargin{1};
        frames=varargin{2};
        ROI=varargin{3};
        ref_ROI=varargin{4};
        thresh=varargin{5};
        flag=varargin{6};
end
%% Import selected frames
mov=extract_frames(filename,frames);

%% View selected frame
% Setup video file for recording
if flag==1
    v=VideoWriter('Sample_004.mp4','MPEG-4');%create video object
    v.FrameRate=3;
    open(v);%open video object
end

for fn=double(frames)
    % Index of frame in mov struct associated with fn
    k=find(mov.abs_frame_index==200);
    
    % Create a figure (rgb)
    f1.f=figure(1); clf(figure(1));
    f1.f.Color='w';
    f1.f.Name=['frame: ',num2str(fn)];
    gap=[0.1 0.04];%Define gap associated with subplot creation
    f1.s1=subtightplot(2,3,[1 4],gap);
    f1.s2=subtightplot(2,3,[2 5],gap);
    f1.s3=subtightplot(2,3,3);
    f1.s4=subtightplot(2,3,6);
    set(f1.s3,'position',[0.7 0.63 0.29 0.27]);
    set(f1.s4,'position',[0.7 0.1 0.29 0.27]);
    set(findall(f1.f,'type','axes'),'nextplot','add','box','on');
    mymap=flipud(parula);%define colormap
    crange=[0 0.05];%c axis range

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define white reference region
    white=mov(k).CData(ref_ROI(1):ref_ROI(2),ref_ROI(3):ref_ROI(4),:);
    white2=rgb_correction(white,white,'simple',250,0);
    [~,~,white_b]=rgb_ratio(white2);
    pd=fitdist(white_b(:),'normal');%fit normal dist to histogram    

    % Extract image frame (rgb)
    rgb_image=mov(k).CData(ROI(1):ROI(2),ROI(3):ROI(4),:);
    rgb_image2=rgb_correction(rgb_image,white,'simple');
    [R_channel,G_channel,B_channel]=rgb_ratio(rgb_image2);
    channel=B_channel;%scale onto grayscale (0 to 255)

    % image correction
    channel2=channel;
    channel3=medfilt2(channel2,[10 10]);
    channel4=rm_lowest_bin(channel3,thresh,0);
    channel5=(channel4-thresh).*(250/255);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtered contour plot with transparent original image overlay
    contour(f1.s1,channel5,'fill','on');
    colormap(f1.s1,mymap);
    f1.p1=imagesc(f1.s1,uint8(rgb_image2));
    f1.p1.AlphaData=1;
    axis(f1.s1,'image');
    colorbar(f1.s1,'southoutside')
    caxis(f1.s1,crange)

    % Filtered contour plot only
    contour(f1.s2,channel5,'fill','on')
    axis(f1.s2,'image');
    colormap(f1.s2,mymap);
    colorbar(f1.s2,'southoutside')
    caxis(f1.s2,crange)

    % Histogram of complment gray scaled image before background subtraction
    histogram(f1.s3,channel2,255,'edgecolor','none',...
        'facecolor','k','facealpha',1);
    f1.p2=plot(f1.s3,[thresh thresh],f1.s3.YLim,'r--');
    title(f1.s3,['Histogram of grayscale complement',char(10),...
        'before background subtraction']);
    xlabel(f1.s3,'Relative Intensity');
    ylabel(f1.s3,'Counts');
    L=legend(f1.p2,['background',char(10),'threshold']);
    L.Color='none';
    f1.s3.XLim=[thresh*0.9 thresh*1.1];

    % Histogram of complement gray scaled image after background subtraction,
    % application of 2D median filter, removal of first bin (background) and
    % spurious pixels
    histogram(f1.s4,channel5,255,'edgecolor','none','facecolor','k',...
        'facecolor','k','facealpha',1);
    title(f1.s4,['Histogram of grayscale',char(10),...
        'image after corrections']);
    xlabel(f1.s4,'Relative Intensity');
    ylabel(f1.s4,'Counts');
    f1.s4.XLim(1)=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % format axes
    % Change tick labels to represent the absolute pixel indices in
    % reference to the original frame.
    [xticklabels,yticklabels]=gen_labels(ROI(3):100:ROI(4),ROI(1):100:ROI(2));
    set(f1.s1,'xtick',0:100:ROI(4)-ROI(3),'xticklabel',xticklabels,...
        'ytick',0:100:ROI(2)-ROI(1),'yticklabel',yticklabels)
    set(f1.s2,'xtick',0:100:ROI(4)-ROI(3),'xticklabel',xticklabels,...
        'ytick',0:100:ROI(2)-ROI(1),'yticklabel',yticklabels)
    set(findall(f1.f,'type','text'),'fontweight','bold');
    linkaxes([f1.s1,f1.s2],'xy');
    
    %Grab the frame from f1.f
    disp(['frame: ',num2str(fn)]);
    
    if flag==1
        F=getframe(f1.f);
        writeVideo(v,F.cdata)    
    end
end

%close video object
if flag==1
    close(v);
end