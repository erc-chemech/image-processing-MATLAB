function recalc_lambda(filename,frames,top_roi,bot_roi,varargin)
%%
% Author: Joshua Yeh
% Date created: 190228
% 
%% DESCRIPTION
% This fcn recalculates the extension ratio, lambda, of a sample that has
% not been marked during a tensile test.
%
%% INPUT VARIABLES
% filename: filename of thevideo to be analyzed
% 
% frames: the frame indices that will be analyzed from the video
% 
% m_file: the mechanical data obtained from an instron
% 
% top_roi: top region of interest that will be tracked marked by a feature
% or dust particle
% 
% bot_roi: bottom region of interest that will be traced marked by a
% feature or dust particle
% 
% NAME PAIR ARGUMENTS: RGB_analysis(...'<fieldname>',<value>)
% 
% 'xlsfile': exported excel file name
% 
% 'm_file': mechanical data obtained from the Instron
% 
% 'rotate': apply rotation to the frame in degrees
% 
%% output variables
% Note: variables in this fcn will be outputted into the caller workspace.
% Also, the time, lambda (extension ratio), and stress will will exported
% to an excel file based on the basename of the video filename
% 
%%
top_roi0=top_roi;
bot_roi0=bot_roi;
%% parse input

narginchk(5,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('xlsfile',[],@(x) ischar(x));
params.addParameter('m_file',[],@(x) isempty(x)|ischar(x));
params.addParameter('rotate',0,@(x) isnumeric(x));
params.parse(varargin{:});

% Extract out values from parsed input
xlsfile=params.Results.xlsfile;
m_file=params.Results.m_file;
rot=params.Results.rotate;

% Turn off hardware acceleration
matlab.video.read.UseHardwareAcceleration('off')
%Create a VideoReader object
Vidobj = VideoReader(filename);
% Turn back on the hardware acceleration
matlab.video.read.UseHardwareAcceleration('on')

% prepare dynamic figures
f1=my_fig(1,{[1 2 1] [1 2 2]});
f1.I=imagesc(f1.s1,uint8(1));
f1.cent=plot(f1.s1,nan,nan,'rx');
f1.pl=plot(f1.s1,nan,nan,'rs');
colormap(bone);
title(f1.s1,'top');
axis(f1.s1,'image');

f2.s1=f1.s2;
f2.I=imagesc(f2.s1,uint8(1));
f2.cent=plot(f2.s1,nan,nan,'rx');
f2.pl=plot(f2.s1,nan,nan,'rs');
colormap(bone);
title(f2.s1,'bottom');
axis(f2.s1,'image');

f3=my_fig(3);
f3.I=imagesc(uint8(1));
f3.top=plot(f3.s1,nan,nan,'rx');
f3.bot=plot(f3.s1,nan,nan,'rx');
axis(f3.s1,'image');

f4=my_fig(4);
xylabels(f4.s1,'Time (s)','\lambda');
f4.p=plot(f4.s1,nan,nan,'o-');
center_axes;

top_dot=nan(numel(frames),2);
bot_dot=nan(numel(frames),2);
dist=nan(numel(frames),1);
time_store=dist;

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
        
        if count2==1%remember first time point
            start=Vidobj.CurrentTime;
        end
        time_store(count2)=Vidobj.CurrentTime-start;%rel. time
        
        frame=mov(1).CData;%orginal frame
        
        % if user defines a rotation or rotation is non-zero, rotate frame
        if rot~=0
            frame=imrotate(frame,rot);
        end

        % Define top portion of the frame
        top0=frame(top_roi(1):top_roi(2),top_roi(3):top_roi(4),:);
        top_offsetx=top_roi(3)-1;
        top_offsety=top_roi(1)-1;
        top1=rgb2gray(top0);
        top2=imcomplement(imbinarize(top1));
        top3=bwareaopen(top2,2);

        % Perform image analysis on the top portion of the frame
        stats_top=regionprops(top3,'Centroid','Area','PixelList');
        areas=[stats_top.Area];
        [~,kk]=max(areas);
        centers_t=stats_top(kk).Centroid;%center coord in top_roi coord system
        centers0_t=[centers_t(:,1)+top_offsetx centers_t(:,2)+top_offsety];
        top_dot(count2,:)=centers0_t;%store the centroid info in array
        PixelList_t=stats_top(kk).PixelList;

        % Define bottom portion of the frame
        bot0=frame(bot_roi(1):bot_roi(2),bot_roi(3):bot_roi(4),:);
        bot_offsetx=bot_roi(3)-1;
        bot_offsety=bot_roi(1)-1;
        bot1=rgb2gray(bot0);
        bot2=imcomplement(imbinarize(bot1));
        bot3=bwareaopen(bot2,2);

        % Perform image analysis on the bottom portion of the frame
        stats_bot=regionprops(bot3,'Centroid','Area','PixelList');
        areas=[stats_bot.Area];
        [~,kk]=max(areas);
        centers_b=stats_bot(kk).Centroid;%center coordinates in top_roi coord system
        centers0_b=[centers_b(:,1)+bot_offsetx centers_b(:,2)+bot_offsety];
        bot_dot(count2,:)=centers0_b;%store the centroid info in array
        PixelList_b=stats_bot(kk).PixelList;

        % Calculate pixel distance between the top and bot dots
        dist(count2)=abs(bot_dot(count2,2)-top_dot(count2,2));
        
        % calculate extension ratio
        lambda=dist./dist(1);

        % update ROI for top and bot
        if count2>1

            % Calculate change in the x direction'
            top_changex=round(top_dot(count2,1)-top_dot(1,1));
            bot_changex=round(bot_dot(count2,1)-bot_dot(1,1));

            % Calculate change in the y direction
            top_changey=round(top_dot(count2,2)-top_dot(1,2));
            bot_changey=round(bot_dot(count2,2)-bot_dot(1,2));


            %update rois
            top_roi(3:4)=top_roi0(3:4)+top_changex;
            top_roi(1:2)=top_roi0(1:2)+top_changey;
            bot_roi(3:4)=bot_roi0(3:4)+bot_changex;
            bot_roi(1:2)=bot_roi0(1:2)+bot_changey;

        end

        % Update plot of analysis
        % figure 1
        set(f1.I,'cdata',top3);
        set(f1.cent,'xdata',centers_t(:,1),'ydata',centers_t(:,2));
        set(f1.pl,'xdata',PixelList_t(:,1),'ydata',PixelList_t(:,2));
        set(f1.s1,'ydir','reverse');

        % figure 2
        set(f2.I,'cdata',bot3);
        set(f2.cent,'xdata',centers_b(:,1),'ydata',centers_b(:,2));
        set(f2.pl,'xdata',PixelList_b(:,1),'ydata',PixelList_b(:,2));
        set(f2.s1,'ydir','reverse');

        % figure 3
        set(f3.I,'cdata',frame);
        set(f3.top,'xdata',centers0_t(:,1),'ydata',centers0_t(:,2));
        set(f3.bot,'xdata',centers0_b(:,1),'ydata',centers0_b(:,2));
        set(f3.s1,'ydir','reverse');

        % figure 4
        set(f4.p,'xdata',time_store,'ydata',lambda);

        drawnow
        count2=count2+1;% update counter
    end
    count1=count1+1;% update counter
end

%%
if ~isempty(m_file)
    % Load the mech data
    [num,txt,raw]=xlsread(m_file);
    %get stress value
    k=find(strcmp(txt,'(MPa)'));
    [~,kc]=ind2sub(size(txt),k);
    m_stress=num(:,kc);

    %get time value
    k=find(strcmp(txt,'(sec)'));
    if isempty(k)
        k=find(strcmp(txt,'(s)'));
    end
    [~,kc]=ind2sub(size(txt),k);
    m_time=num(:,kc);

    %interpolate stress data on vid time scale
    stress_q=spline(m_time,m_stress,time_store);
else
    stress_q=nan(size(time_store));
end

%export results in excel spreadsheet
data=[time_store,lambda,stress_q,frames'];
data=mat2cell(data,ones(size(data,1),1),[1 1 1 1]);
headers={'Time (s)','Lambda','Stress (MPa)','frame #'};

if isempty(xlsfile)
    [filepath,name,ext] = fileparts(filename);
    name=[name,'_recalc_lam'];
    xlswrite(name,[headers;data]);
else
    xlswrite(xlsfile,[headers;data]);
end

% export fcn workspace to caller workspace
W=who;
out_var(W{:});

disp(['fini: ',datestr(clock)]);

function out_var(varargin)
% This function output the function variable space to the base workspace
for dum=1:numel(varargin)
    assignin('base',varargin{dum},evalin('caller',varargin{dum}));
end