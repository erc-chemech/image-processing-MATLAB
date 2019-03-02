function RGB_analysis(filename,frames,varargin)
% Author: Joshua Yeh
% Date created: 19-02-27
% 
%% DESCRIPTION
% The script performs rgb analysis of a uniaxial sample containing
% spiropyran. The sample must have black lines drawn actoss the sample.
% This script calcualtes the change in chromaticity for the red, green, and
% blue channels.
% 
%% INPUT VARIABLES
% filename: filename of the tensile test video that will be analyzed
% 
% frames: the frame indices to perform the analysis
% 
% NAME PAIR ARGUMENTS: RGB_analysis(...'<fieldname>',<value>)
% 
% 'ROI': region of interest within the original frame
% 
% 'white_ROI': region of interest within the cropped ROI representing the
% reference white region
% 
% 'black_thresh': threshold value that identifies the dark areas
% 
% 'h_divider': horizontal line dividing the sample (rel. to cropped ROI)
% 
% 'v_divider': vertical line dividing the sample (del. to cropped ROI)
% 
% 'mech_file': file to the mechanical data (.xlsx)
% 
% 'export_filename': the filename for the exported results (default is
% based on the base name of the video file)
% 
% 'rotate': rotate counterclockwise direction the frame in degrees
% 
%% OUTPUT
% Note that all the variables in this fcn will be placed in the caller
% workspace.

%% PARSE INPUT VARIABLES

narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('ROI',[550 1200 430 550],@(x) isnumeric(x));
params.addParameter('white_ROI',[450 550 10 25],@(x) isnumeric(x));
params.addParameter('black_thresh',40,@(x) isnumeric(x));
params.addParameter('h_divider',500,@(x) isnumeric(x));
params.addParameter('v_divider',65,@(x) isnumeric(x));
params.addParameter('mech_file',[],@(x) ischar(x));
params.addParameter('export_filename',[],@(x) ischar(x));
params.addParameter('rotate',[],@(x) isempty(x)|isnumeric(x));
params.parse(varargin{:});

% Extract out values from parsed input
ROI=params.Results.ROI;
white_ROI=params.Results.white_ROI;
black_thresh=params.Results.black_thresh;
h_divider=params.Results.h_divider;
v_divider=params.Results.v_divider;
m_file=params.Results.mech_file;
export_filename=params.Results.export_filename;
rot=params.Results.rotate;

% Preallocate arrays to speed algorithm
abs_frame_index=nan(numel(frames),1);% absolute frame number
R_mean=nan(2e4,1);%R-ratio mean
G_mean=nan(2e4,1);%G-ratio mean
B_mean=nan(2e4,1);%B-ratio mean
R_std=nan(2e4,1);%R-ratio stdev
G_std=nan(2e4,1);%R-ratio stdev
B_std=nan(2e4,1);%R-ratio stdev
time_store=nan(2e4,1);% time associated with the frame
lambda=nan(2e4,1);%extension ratio
dist=nan(2e4,1);%distance between the black lines

if ~isempty(m_file)
    [num,txt,raw]=xlsread(m_file);
    
    %get stress value
    k=find(strcmp(txt,'(MPa)'));
    [~,kc]=ind2sub(size(txt),k);
    m_stress=num(:,kc);
    
    %get time value
    k=find(strcmp(txt,'(sec)'));
    [~,kc]=ind2sub(size(txt),k);
    m_time=num(:,kc);
end
        
% Figure creation
f1=my_fig(1,{[3 2 1 5] [3 2 2] [3 2 6] [3 2 4]},'fontsize',14);

% Figure and axes formatting
set(f1.s1,'fontsize',6,'ydir','reverse','box','off');
axis(f1.s1,'image');
set(findall(f1.f,'type','axes'),'nextplot','add');

% Create inital plot handles
f1.image=imagesc(f1.s1,uint64(ones(2)));
f1.horizontal=plot(f1.s1,nan,nan,'k--');
f1.vertical=plot(f1.s1,nan,nan,'k--');
f1.ROI=plot(f1.s1,nan,nan,'r-');
f1.white=plot(f1.s1,nan,nan,'r--');
f1.top_bar=plot(f1.s1,nan,nan,'r-');
f1.bot_bar=plot(f1.s1,nan,nan,'r-');

f1.mR=plot(f1.s3,nan,nan,'r.');
f1.mG=plot(f1.s3,nan,nan,'g.');
f1.mB=plot(f1.s3,nan,nan,'b.');
xylabels(f1.s3,'nominal stress (MPa)','\DeltaRGB');

f1.R_mean=plot(f1.s2,time_store,R_mean,'r.','markersize',14);
f1.G_mean=plot(f1.s2,time_store,R_mean,'g.','markersize',14);
f1.B_mean=plot(f1.s2,time_store,R_mean,'b.','markersize',14);
xylabels(f1.s2,'time (s)','RGB ratio');

f1.delR=plot(f1.s4,nan,nan,'r.','markersize',14);
f1.delG=plot(f1.s4,nan,nan,'g.','markersize',14);
f1.delB=plot(f1.s4,nan,nan,'b.','markersize',14);
xylabels(f1.s4,'time (s)','\DeltaRGB ratio');

% Plot white reference area
set(f1.white,'xdata',white_ROI([3 3 4 4 3]),...
    'ydata',white_ROI([1 2 2 1 1]));

%%
% Video file is too large to import all of the frames into memory. Thus,
% individual frames will be imported one at a time through the for loop.

% Turn off hardware acceleration
matlab.video.read.UseHardwareAcceleration('off')
%Create a VideoReader object
Vidobj = VideoReader(filename);
% Turn back on the hardware acceleration
matlab.video.read.UseHardwareAcceleration('on')

% Import frames from video into mov struct variable
mov = struct('CData',zeros(Vidobj.Height,Vidobj.Width,3,1,'uint8'),...
    'abs_frame_index',[],'CurrentTime',[]);
mov(1).CData=[];%preallocate structure array

flag1=1;% This flag controls the condition in which to continue while loop
count1=1;% Counter counting the iteration for each while loop
count2=1;% Counter counting the iteration # in the frame array being looped
check=0;
while hasFrame(Vidobj)&&count1<=max(frames)
    
    % Import current frame
    mov(1).CData=readFrame(Vidobj);%read and store current frame
    mov(1).abs_frame_index=count1;%abs frame index
    mov(1).CurrentTime=Vidobj.CurrentTime;%abs frame time
    
    % Prevent processing of duplicate frames
    if count1>1
        if isequal(prev,mov(1).CData)&&count2>1
            check=0;
            disp('Skipped duplicate frame');
            
        else
            check=1;
        end
    end
    prev=mov(1).CData;
    
    %check to see whether or not to process frame
    if sum(frames==count1)==1&&check==1
        disp(['Current frame time: ',num2str(Vidobj.CurrentTime),...
            ' ',num2str(count1)]);
        abs_frame_index(count2)=count1;% abs store frame index
        
        if count2==1%remember first time point
            start=Vidobj.CurrentTime;
        end
        time_store(count2)=Vidobj.CurrentTime-start;%rel. time
        
        % interpolate stresses if mech file is provided
        if ~isempty(m_file)
            stress_q=spline(m_time,m_stress,time_store);
        end
        
        % Crop out ROI
        frame=mov(1).CData(ROI(1):ROI(2),ROI(3):ROI(4),:);
        if ~isempty(rot)%check to see if frame needs to be rotate
            frame=imrotate(frame,rot);%rotate the frame
            frame=frame(20:end-20,40:end-40,:);%remove edges caused by rotation
        end
        
        % Define white reference area
        white=frame(white_ROI(1):white_ROI(2),...
            white_ROI(3):white_ROI(4),:);
        
        % Apply color correction
        frame_corr=rgb_correction(frame,white,'simple',200);
        
        b_channel=frame_corr(:,:,3);%blue channel
        b_channel1=medfilt2(b_channel,[5 5]);%apply 2d median filter fcn

        ii1=find(b_channel1(:)<=black_thresh);
        [r1,~]=ind2sub(size(b_channel1),ii1);%extract row
        [r1,~]=sort(r1);%sort the rows in asending order

        % Find top boundary
        ii3=find(r1<h_divider,1,'last');
        top_black=r1(ii3)+1;

        % Find bottom boundary
        ii3=find(r1>h_divider,1,'first');
        bottom_black=r1(ii3)-1;
        
        % Calculate distance between top and bottom
        dist(count2)=abs(top_black-bottom_black);
        
        % Calculate the extension ratio, lambda
        lambda(count2)=dist(count2)./dist(1);

        % Define subimage defined by top and bottom boundariers that will be
        % used to find the left and right boundaries
        b_subimage=b_channel1(top_black-10:bottom_black+10,:);

        % Find indices corresponding to pixels falling wihtin threshold of the
        % subimage
        ii1=find(b_subimage(:)<=black_thresh);

        [~,c1]=ind2sub(size(b_subimage),ii1);% extract row
        [c1,~]=sort(c1);

        % Find left boundary
        ii3=find(c1<=v_divider,1,'first');
        left_black=c1(ii3)+5;

        % Find right boundary
        ii3=find(c1>=v_divider,1,'last');
        right_black=c1(ii3)-5;

        %Define subimage that will be used to calculate statistics
        subimage=frame_corr(top_black+15:bottom_black-15,...
            left_black:right_black,:);
        
        subimage0=frame(top_black+15:bottom_black-15,...
            left_black:right_black,:);
        
        % Determine the RGB ratio of the subimage and calculate stats
        [R_ratio,G_ratio,B_ratio]=rgb_ratio(subimage);
        R_mean(count2)=nanmean(R_ratio(:));
        G_mean(count2)=nanmean(G_ratio(:));
        B_mean(count2)=nanmean(B_ratio(:));
        R_std(count2)=nanstd(R_ratio(:));
        G_std(count2)=nanstd(G_ratio(:));
        B_std(count2)=nanstd(B_ratio(:));
        
        % Based on the top boundary line of the subimage, recalculate
        % appropriate h_divider(n)
        if count2>1
            h_divider=bottom_black-round((bottom_black-top_black)/2);
            v_divider=right_black-round((right_black-left_black)/2);

            % Calculate the Delta RGB ratio
            RCC=R_mean-nanmean(R_mean(1:5));
            GCC=G_mean-nanmean(G_mean(1:5));
            BCC=B_mean-nanmean(B_mean(1:5));

            % Plot delta_RGB values
            set(f1.delR,'xdata',time_store,...
                'ydata',RCC);
            set(f1.delG,'xdata',time_store,...
                'ydata',GCC);
            set(f1.delB,'xdata',time_store,...
                'ydata',BCC);
            
            % Plot chromaticity vs mechanical response
            if ~isempty(m_file)
                set(f1.mR,'xdata',stress_q,'ydata',RCC);
                set(f1.mG,'xdata',stress_q,'ydata',GCC);
                set(f1.mB,'xdata',stress_q,'ydata',BCC);
            end
        end
        
        % Display the color corrected image
        set(f1.image,'cdata',uint8(rgb_correction(frame,white,'simple',255)));

        % Plot horizontal divider    
        set(f1.horizontal,'XData',[1 size(frame_corr,2)],...
            'YData',ones(1,2).*(h_divider));

        % Plot vertical divider    
        set(f1.vertical,'xdata',ones(1,2).*v_divider,...
            'ydata',[1 size(frame_corr,1)]);
        
        % Plot ROI
        set(f1.ROI,'xdata',...
            [left_black right_black right_black left_black left_black],...
            'ydata',...
            [bottom_black-10 bottom_black-10 top_black+10 top_black+10,...
            bottom_black-10]);
        set(f1.top_bar,'xdata',[left_black right_black],...
            'ydata',[top_black top_black]);
        set(f1.bot_bar,'xdata',[left_black right_black],...
            'ydata',[bottom_black bottom_black]);

        % Plot mean of subimage area
        set(f1.R_mean,'YData',R_mean,'xdata',time_store);
        set(f1.G_mean,'YData',G_mean,'xdata',time_store);
        set(f1.B_mean,'YData',B_mean,'xdata',time_store);
        
        
        drawnow;
        count2=count2+1;% update counter
    end
    count1=count1+1;% update counter
end

linkaxes([f1.s2 f1.s4],'x');

% get rid of extraneous nans at the end of the array
R_mean=R_mean(~isnan(time_store));
G_mean=G_mean(~isnan(time_store));
B_mean=B_mean(~isnan(time_store));
R_std=R_std(~isnan(time_store));
G_std=G_std(~isnan(time_store));
B_std=B_std(~isnan(time_store));
RCC=RCC(~isnan(time_store));
GCC=GCC(~isnan(time_store));
BCC=BCC(~isnan(time_store));
lambda=lambda(~isnan(time_store));
time_store=time_store(~isnan(time_store));

% export fcn workspace to caller workspace
W=who;
out_var(W{:});

if isempty(export_filename)
    [filepath,name,ext] = fileparts(filename);
    export_filename=[name,'_results.mat'];
end
if isempty(m_file)
    save(export_filename,'time_store','R_mean','G_mean','B_mean',...
        'RCC','GCC','BCC','R_std','G_std','B_std','lambda','white',...
        'subimage','subimage0','varargin');
else
    save(export_filename,'time_store','R_mean','G_mean','B_mean',...
        'RCC','GCC','BCC','R_std','G_std','B_std','lambda','white',...
        'subimage','subimage0','varargin','stress_q','m_stress','m_time',...
        'raw');
end

disp(['Finished on: ',datestr(clock)])
disp(['results saved to: ',export_filename]);

function out_var(varargin)
% This function output the function variable space to the base workspace
for dum=1:numel(varargin)
    assignin('base',varargin{dum},evalin('caller',varargin{dum}));
end
