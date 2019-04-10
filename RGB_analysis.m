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
% 'ROI': region of interest within the original frame (prior to frame
% rotation)
% 
% 'white_ROI': region of interest with respect to the original (uncropped)
% frame prior to rotation
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
% 'border_offset': border offset when identifying the sample region
% 
% 'crop_border': crop the border after rotation (will only apply if the
% value for 'rotate' ~=0) (pxs)
% 
% 'time_offset': define a time offset (default is 0)
% 
% 'CC_ref': define ref for calc. chromatic changes (default is [0 0
% 0], corresponding to R, G, B, respectively). Note that the value array
% must have three elements.
% 
%% OUTPUT
% Note that all the variables in this fcn will be placed in the caller
% workspace.

%% PARSE INPUT VARIABLES

narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('ROI',[],@(x) isnumeric(x));
params.addParameter('white_ROI',[450 550 10 25],@(x) isnumeric(x));
params.addParameter('black_thresh',150,@(x) isnumeric(x));
params.addParameter('h_divider',500,@(x) isnumeric(x));
params.addParameter('v_divider',65,@(x) isnumeric(x));
params.addParameter('mech_file',[],@(x) ischar(x));
params.addParameter('export_filename',[],@(x) ischar(x));
params.addParameter('rotate',0,@(x) isnumeric(x));
params.addParameter('border_offset',15,@(x) isnumeric(x));
params.addParameter('crop_border',50,@(x) isnumeric(x));
params.addParameter('time_offset',0,@(x) isnumeric(x));
params.addParameter('CC_ref',[0 0 0],@(x) isnumeric(x)&numel(x)==3);
params.parse(varargin{:});

% Extract out values from parsed input
ROI=params.Results.ROI;
white_ROI=params.Results.white_ROI;
black_thresh=params.Results.black_thresh;
h_divider=params.Results.h_divider;
v_divider=params.Results.v_divider;
m_file=params.Results.mech_file;
export_filename=params.Results.export_filename;
theta=params.Results.rotate;
border_offset=params.Results.border_offset;
crop_border=params.Results.crop_border;
t_offset=params.Results.time_offset;
CC_ref=params.Results.CC_ref;

% Preallocate arrays to speed algorithm
abs_frame_index=nan(numel(frames),1);% absolute frame number
R_mean=nan(2e4,1);%R-ratio mean
G_mean=nan(2e4,1);%G-ratio mean
B_mean=nan(2e4,1);%B-ratio mean
R_std=nan(2e4,1);%R-ratio stdev
G_std=nan(2e4,1);%R-ratio stdev
B_std=nan(2e4,1);%R-ratio stdev
time_store=nan(2e4,1);% time associated with the video frame
lambda=nan(2e4,1);%extension ratio
dist=nan(2e4,1);%distance between the black lines

% extract mechanical data
disp('loading mechanical data');
if ~isempty(m_file)
    [~,txt,raw]=xlsread(m_file);
    
    %get stress value
    k=find(strcmp(txt,'(MPa)'));
    [kr,kc]=ind2sub(size(txt),k);
    m_stress=[raw{kr+1:end,kc}]';
    
    %get time value
    k=find(strcmp(txt,'(sec)'));
    if isempty(k)
        k=find(strcmp(txt,'(s)'));
    end
    [kr,kc]=ind2sub(size(txt),k);
    m_time=[raw{kr+1:end,kc}]';
end
        
% Figure creation
f1=my_fig(1,{[3 2 1 5] [3 2 2] [3 2 6] [3 2 4]},'fontsize',14);

% Figure and axes formatting
set(f1.s1,'fontsize',6,'ydir','reverse','box','off');
axis(f1.s1,'image');
set(findall(f1.f,'type','axes'),'nextplot','add');

% Create dynamic plot handles

% axes associated with the frame
f1.image=imagesc(f1.s1,uint64(ones(2)));
f1.horizontal=plot(f1.s1,nan,nan,'k--');
f1.vertical=plot(f1.s1,nan,nan,'k--');
f1.ROI=plot(f1.s1,nan,nan,'r-');
f1.white=plot(f1.s1,nan,nan,'r--');
f1.top_bar=plot(f1.s1,nan,nan,'r-');
f1.bot_bar=plot(f1.s1,nan,nan,'r-');

% chromatic values vs time
f1.R_mean=plot(f1.s2,time_store,R_mean,'r.','markersize',14);
f1.G_mean=plot(f1.s2,time_store,R_mean,'g.','markersize',14);
f1.B_mean=plot(f1.s2,time_store,R_mean,'b.','markersize',14);
xylabels(f1.s2,'time (s)','RGB ratio');

% chromatic change vs stress
f1.mR=plot(f1.s3,nan,nan,'r.');
f1.mG=plot(f1.s3,nan,nan,'g.');
f1.mB=plot(f1.s3,nan,nan,'b.');
xylabels(f1.s3,'nominal stress (MPa)','\DeltaRGB');

% chromatic change vs time
f1.delR=plot(f1.s4,nan,nan,'r.','markersize',14);
f1.delG=plot(f1.s4,nan,nan,'g.','markersize',14);
f1.delB=plot(f1.s4,nan,nan,'b.','markersize',14);
xylabels(f1.s4,'time (s)','\DeltaRGB ratio');

%%
% Video file is too large to import all of the frames into memory. Thus,
% individual frames will be imported one at a time through the for loop.

disp('creating video object');
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

count1=1;% Counter counting the iteration for each while loop
count2=1;% Counter counting the iteration # in 'frame' array being looped
check=0;% a flag for checking if frame is new (0 for dupli. 1 for new)
while hasFrame(Vidobj)&&count1<=max(frames)
    
    % Import current frame
    mov(1).CData=readFrame(Vidobj);%read and store current frame
    mov(1).abs_frame_index=count1;%abs frame index
    mov(1).CurrentTime=Vidobj.CurrentTime;%abs frame time
    
    % Prevent processing of duplicate frames
    if count1>1
        if isequal(prev,mov(1).CData)&&count2>1
            check=0;% mark as duplicate frame
            disp('Skipped duplicate frame');
            
        else
            check=1;% mark as a new frame
        end
    elseif count1==1
        check=1;% for the 1st frame, toggle check as 1 (new frame)
    end
    prev=mov(1).CData;
    
    % mark the white ROI (only need to do this once)
    if count1==1
        
        % Create logical array w/ the size of the raw frame. This will be
        % used to track the white_ROI location during frame transformation
        % (after cropping and rotation).
        marked=false(size(mov(1).CData(:,:,1)));
        marked(white_ROI(1),white_ROI(3))=true;%mark upper left pt
        marked(white_ROI(1),white_ROI(4))=true;%mark upper right pt
        marked(white_ROI(2),white_ROI(4))=true;%mark lower right pt
        marked(white_ROI(2),white_ROI(3))=true;%mark lower left pt
        
        % Crop out ROI from raw image if user specifies ROI
        if ~isempty(ROI)%if ROI is specified
            marked=marked(ROI(1):ROI(2),ROI(3):ROI(4),:);
        end
        
        % check to see if frame needs to be rotate
        if theta~=0
            marked=imrotate(marked,theta);%rotate logical array
            
            %remove edges due to rotation
            marked=marked(crop_border+1:end-crop_border,...
                crop_border+1:end-crop_border,:);
        end
        
        % Determine the white_ROI coordinates after cropping & rotation
        mark_idx=find(marked==true);
        [mark_y,mark_x]=ind2sub(size(marked),mark_idx);
        mark_xy=[mark_x(:),mark_y(:)];%xy coordinates of marked pts
        
        % Need to perform knnsearch (nearest neighbor search) so that the
        % order of the xy-coordinates will result in the plotting of a box.
        mark_xy_idx=knnsearch(mark_xy,mark_xy(1,:),'k',4);
        mark_xy2=mark_xy(mark_xy_idx([1;2;4;3]),:);%reorder pts
        mark_xy2=cat(1,mark_xy2,mark_xy2(1,:));%close the box
        
        % show the white_ROI
        set(f1.white,'xdata',mark_xy2(:,1),'ydata',mark_xy2(:,2));
    end
    
    %check to see whether or not to process frame
    if sum(frames==count1)==1&&check==1
        disp(['Current frame time: ',num2str(Vidobj.CurrentTime),...
            ' ',num2str(count1)]);
        abs_frame_index(count2)=count1;% abs store frame index
        
        if count2==1%remember first time point
            start=Vidobj.CurrentTime;
        end
        time_store(count2)=Vidobj.CurrentTime-start+t_offset;%rel. time
        
        % interpolate stresses if mech file is provided
        if ~isempty(m_file)
            stress_q=spline(m_time,m_stress,time_store);
        end
        
        frame=mov(1).CData;%extract current raw frame
        
        % Define white reference area from raw image
        white=mov(1).CData(white_ROI(1):white_ROI(2),...
            white_ROI(3):white_ROI(4),:);        
        
        % Crop out ROI from raw image if user specifies ROI
        if ~isempty(ROI)%if ROI is specified
            frame=frame(ROI(1):ROI(2),ROI(3):ROI(4),:);
        end
        
        % check to see if frame needs to be rotate
        if theta~=0
            frame=imrotate(frame,theta);%rotate the frame
            
            %remove edges due to rotation
            frame=frame(crop_border+1:end-crop_border,...
                crop_border+1:end-crop_border,:);
        end
                  
        % Apply color correction
        frame_corr=rgb_correction(frame,white,'simple',200,'flag',0);
        
        % define channel1, which represents a subimage within 'frame'
        % defined by 5 pixels to the left and to the right of the vertical
        % dividing line. The subimage is then summed in the 3rd dimension.
        % In other words, it is an intensity array.
        channel1=sum(frame_corr(:,v_divider-5:v_divider+5,:),3);
        channel1b=nanmean(channel1,2);%take the mean in the col direction
        
        % find indices that are equal or below the black_thresh value
        ii1=find(channel1b(:)<=black_thresh);

        % Find top boundary
        ii3=find(ii1<h_divider,1,'last');
        top_black=ii1(ii3)+1;

        % Find bottom boundary
        ii3=find(ii1>h_divider,1,'first');
        bottom_black=ii1(ii3)-1;
        
        % Calculate distance between top and bottom
        try
            dist(count2)=abs(top_black-bottom_black);
        catch
            % show histogram of intensity subimage
            figure; histogram(channel1b(:));
            xylabels(gca,'intensity (a.u.)','counts');
            center_axes;
            
            % show subimage
            figure; imagesc(frame); axis image; center_axes;
            
            error(['Unable to continue. Value for black_thresh is ',...
                'likely badly conditioned. A histogram is plotted' ,...
                'to visualize appropriate black_thresh value. It ',...
                'is possible that the h_divider and/or v_divider ',...
                'is also badly conditioned. An image of the frame ',...
                'is shown to visualize approprate value(s) ',...
                'for h_divider and/or v_divider.']);
        end
        
        % Calculate the extension ratio, lambda
        lambda(count2)=dist(count2)./dist(1);

        % Define subimage defined by top and bottom boundaries that will be
        % used to find the left and right boundaries
        
        % Calculate the intensity of the subimage
        b_subimage=sum(frame_corr(top_black-10:bottom_black+10,:,:),3);

        % Find indices corresponding to pixels falling within threshold of the
        % subimage
        ii1=find(b_subimage(:)<=black_thresh);

        [~,c1]=ind2sub(size(b_subimage),ii1);% extract row
        [c1,~]=sort(c1);

        % Find left boundary
        ii3=find(c1<=v_divider,1,'first');
        left_black=c1(ii3)+border_offset;

        % Find right boundary
        ii3=find(c1>=v_divider,1,'last');
        right_black=c1(ii3)-border_offset;

        %Define subimage that will be used to calculate statistics
        subimage=...
            frame_corr(top_black+border_offset:bottom_black-border_offset,...
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
        % appropriate h_divider and v_divider
        if count2>1
            h_divider=bottom_black-round((bottom_black-top_black)/2);
            
            % estimate what the vd_test should be
            v_divider=right_black-round((right_black-left_black)/2);

            % Calculate the Delta RGB ratio
            if sum(CC_ref)==0
                RCC=R_mean-nanmean(R_mean(1:5));
                GCC=G_mean-nanmean(G_mean(1:5));
                BCC=B_mean-nanmean(B_mean(1:5));
            else
                RCC=R_mean-CC_ref(1);
                GCC=G_mean-CC_ref(2);
                BCC=B_mean-CC_ref(3);
            end
            
            % total chromatic change (Euclidean norm)
            TCC=sqrt(RCC.^2+GCC.^2+BCC.^2);

            % Plot delta_RGB values
            set(f1.delR,'xdata',time_store,'ydata',RCC);
            set(f1.delG,'xdata',time_store,'ydata',GCC);
            set(f1.delB,'xdata',time_store,'ydata',BCC);
            
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
        
        % Plot sample ROI
        set(f1.ROI,'xdata',...
            [left_black right_black right_black left_black left_black],...
            'ydata',...
            [bottom_black-border_offset bottom_black-border_offset,...
            top_black+border_offset top_black+border_offset,...
            bottom_black-border_offset]);
        
        % Plot sample top line
        set(f1.top_bar,'xdata',[left_black right_black],...
            'ydata',[top_black top_black]);
        
        % Plot sample bottom line
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
if ~isempty(m_file)
    stress_q=stress_q(~isnan(time_store));
end
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
        'subimage','varargin','TCC');
else
    save(export_filename,'time_store','R_mean','G_mean','B_mean',...
        'RCC','GCC','BCC','R_std','G_std','B_std','lambda','white',...
        'subimage','varargin','stress_q','m_stress','m_time',...
        'raw','TCC');
end

disp(['Finished on: ',datestr(clock)])
disp(['results saved to: ',export_filename]);

%% USEFUL FCNs
function out_var(varargin)
% This function output the function variable space to the base workspace
for dum=1:numel(varargin)
    assignin('base',varargin{dum},evalin('caller',varargin{dum}));
end
