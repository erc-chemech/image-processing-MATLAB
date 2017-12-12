function mov=extract_frames(filename,frames,varargin)
% This functions serves to extract a range of frames from a video file. The
% frames are stored in a structure variable.
% Author: Joshua Yeh
% Date created: 2017/10/28
%
%% SYNTAX
% mov=extract_frames(filename,frames)
% mov=extract_frames(filename,'time interval',time_interval)
% mov=extract_frames(filename,frames,'ROI',ROI)
% mov=extract_frames(filename,'time interval',time_interval,'ROI',ROI)
%% INPUT VARIABLES
% 
% mov=extract_frames(filename,frames)
% 
% filename: name of video file (RGB file)
% 
% frames: frame indices that will be extracted from the video file (single
% row or col array containing integers or of type uint64).
% 
%   -If user inputs the string 'all' for the frame_range value, the
%   function will import all of the frames from the video file.
% 
% 
% 
% mov=extract_frames(filename,'time interval',time_interval)
% 
% filename: name of video file (RGB file)
% 
% time_interval: a 3 element double array describing the time interval for
% extracting frames. 
% 
%       time_interval(1): The first element is the start time
% 
%       time_interval(2): The second element is the end time.
% 
%       time_interval(3): The third element describes the frame period in
%       which to extract a frame. For example:
% 
%               if time_interval(3)=1, then collect every frame within the
%               time interval
%               if time_interval(3)=5, then collect every 5 frames within
%               the time interval
% 
% 
% 
% mov=extract_frames(filename,frames,'ROI',ROI)
% 
% filename: name of video file (RGB file)
% 
% frames: frame indices that will be extracted from the video file (single
% row or col array containing integers or of type uint64).
% 
%   -If user inputs the string 'all' for the frame_range value, the
%   function will import all of the frames from the video file.
% 
% ROI: a 4 element double array representing a region of interest (ROI)
% 
%       ROI(1): left boundary
% 
%       ROI(2): right boundary
% 
%       ROI(3): bottom boundary
% 
%       ROI(4): top boundary
% 
%% OUTPUT VARIABLES
% mov: a structure variable containing the extracted frames
%   fieldnames:
%       CData: frame data
%       abs_frame_index: absolute frame index (in reference to video file)
%       CurrentTime: the time stamp associated with the frame

% Turn off hardware acceleration to prevent VideoReader from crashing (an
% apparant graphics card issue that crashes MATLAB)

% Specify default variable values

flag_ROI=0;

switch nargin
    case 1 %not enough inputs
        disp('Frame import aborted!');
        error('Not enough arguments!');
    case 2 %user specified frame variable
        if strcmp(frames,'all')==1&&ischar(frames)==1
            disp('This fcn will attempt to import all frames.');
        elseif (strcmp(frames,'all')==0||strcmp(frames,'time interval')==0)&&...
            ischar(frames)==1
        disp('Frame import aborted!');
        error(['Did not recognize char input for input var ''frames''. ',...
            '\n''frames'' appears to be a %s'],class(frames));
        end
    case 3 %user specified time interval to import frames
        if strcmp(frames,'time interval')==1&&ischar(frames)==1
               %check that the user properly input time interval
                if isa(varargin{1},'double')&&length(varargin{1})==3
                    %time interval
                    time_interval=[sort(varargin{1}(1:2)) varargin{1}(3)];
                    disp(['Extracting frames from time ',num2str(time_interval(1)),...
                        ' to ',num2str(time_interval(2)),'.']);
                    disp(['Data is extracted every ',num2str(time_interval(3)),...
                        ' frame(s).']);
                else
                    disp('Frame import aborted!');
                    error(['''time_interval'' is not recognized to be a',...
                        '2 element double array']);
                end  
        end
    case 4 %user specified frame range and ROI
        if strcmp(varargin{1},'ROI')==1
            ROI=varargin{2};
            % Check to see if the elements have been define properly
            if ROI(2)<ROI(1)||ROI(4)<ROI(3)
                disp('Frame import aborted!');
                error(['ROI varargin not defined properly!',...
                    newline,'Make sure that ROI(2)>ROI(1) and',...
                    'ROI(4)>ROI(3)']);
            end
            flag_ROI=1;
            disp('ROI detected!');
        else
            disp('Frame import aborted!');
            error('Unknown varargin (fcn syntax incorrect).');
        end
    case 5 %user specified time itnerval and ROI
        if strcmp(frames,'time interval')==1&&ischar(frames)==1
               %check that the user properly input time interval
                if isa(varargin{1},'double')&&length(varargin{1})==3
                    %time interval
                    time_interval=[sort(varargin{1}(1:2)) varargin{1}(3)];
                    disp(['Extracting frames from time ',num2str(time_interval(1)),...
                        ' to ',num2str(time_interval(2)),'.']);
                    disp(['Data is extracted every ',num2str(time_interval(3)),...
                        ' frame(s).']);
                else
                    disp('Frame import aborted!');
                    error(['''time_interval'' is not recognized to be a',...
                        '2 element double array']);
                end  
        end
        if strcmp(varargin{2},'ROI')==1
            ROI=varargin{3};
            % Check to see if the elements have been define properly
            if ROI(2)<ROI(1)||ROI(4)<ROI(3)
                disp('Frame import aborted!');
                error(['ROI varargin not defined properly!',...
                    newline,'Make sure that ROI(2)>ROI(1) and',...
                    'ROI(4)>ROI(3)']);
            end
            flag_ROI=1;
            disp('ROI detected!');
        else
            disp('Frame import aborted!');
            error('Unknown varargin (fcn syntax incorrect).');
        end
end

% Check to see if the frames array contains integers
if isinteger(frames)==0&&ischar(frames)==0
    disp('Variable input for frames is not of type uint64.');
    disp('Attempting to convert ''frames'' into uint64.');
    try
        frames=uint64(frames);
    catch
        disp('Unable to convert to uint8. Frame import aborted.');
        return
    end
end

% Turn off hardware acceleration
matlab.video.read.UseHardwareAcceleration('off')

%Create a VideoReader object
Vidobj = VideoReader(filename);

% Turn back on the hardware acceleration
matlab.video.read.UseHardwareAcceleration('on')

% Import frames from video into mov struct variable
mov = struct('CData',zeros(Vidobj.Height,Vidobj.Width,3,1,'uint8'),...
    'abs_frame_index',[],'CurrentTime',[]);
mov(5000).CData=[];%preallocate structure array

n=1;%counter absolute frame index for the video file
k=1;%counter index for mov structure
kk=1;%counter used when user designates time interval (to det. frame extraction)

disp('Importing frames (this can take awhile).');
while hasFrame(Vidobj)
    % If the frame falls within the designated frame/time range, store the
    % frame in the mov structure variable.
    
    %if user specifies range of frames to import
    if ischar(frames)==0
        if ismember(n,frames)==1
            flag=1;
        elseif n<max(frames) % Otherwise, do not store the frame.
            flag=0;
        elseif n>max(frames)            
            break
        end
        
    %If user specifies to import frames within a specified time interval
    elseif ischar(frames)==1&&strcmp(frames,'time interval')
        if Vidobj.CurrentTime>=time_interval(1)&&...
                Vidobj.CurrentTime<time_interval(2)
            if mod(kk,time_interval(3))==1||time_interval(3)==1
                flag=1;
            else
                flag=0;
            end
            kk=kk+1;
        elseif Vidobj.CurrentTime<time_interval(1)
            flag=0;            
        else            
            break;
        end
        
    % If user specifies to import all frames.
    elseif ischar(frames)==1&&strcmp(frames,'all')
            flag=1;
    end
    
    if flag==1%if flag is 1, extract frame
        mov(k).CData=readFrame(Vidobj);%store the frame
        
        % Check if the user specified a ROI
        if flag_ROI==1        
            mov(k).CData=mov(k).CData(ROI(1):ROI(2),ROI(3):ROI(4),:);
        end
        
        mov(k).abs_frame_index=n;%store the absolute frame index
        mov(k).CurrentTime=Vidobj.CurrentTime;%store associated timepoint
        
        if k>1%this is to prevent duplicate frames being saved
            if isequal(mov(k).CData,mov(k-1).CData)==0
                if mod(k,100)==0%update user every 100 frames are read
                    disp([num2str(k),' frames imported.']);
                    disp(['Current frame time: ',num2str(Vidobj.CurrentTime)]);
                end
                k=k+1;
            else                
            end
        else
            k=k+1;
        end        
    else%otherwise advance to the next frame
        readFrame(Vidobj);
        if mod(n,100)==0%update user every 100s on movie timeline
            disp(['Played through ',num2str(n),' frames']);
            disp(['Current frame time: ',num2str(Vidobj.CurrentTime)]);
        end
    end
    
    % Warn the user if the number of extracted frames exceed 500
    if mod(k,100)==0&&k>500
        [user,~] = memory;%get memory
        if user.MemAvailableAllArrays/1e6<3000
            % display warning and memory stats
            warning([newline,'MemAvailableAllArrays: ',...
                num2str(user.MemAvailableAllArrays./1e6),'MB',newline,...
                'MemUsedMATLAB: ',num2str(user.MemUsedMATLAB./1e6),'MB',newline,...
                'Number of frames extracted is over 500. '...
                'MATLAB can crash if memory usage is too high!',newline,...
                'Use Crtl+c in the command window to terminate!']);
            pause(2);
        end
    end
    n=n+1;
end

% Remove empty unused preallocated entries
mov=mov(1:k-1);

disp('Import successful');
disp(['Imported ',num2str(k-1),' frame(s)']);
[user,~] = memory;
disp([newline,'MemAvailableAllArrays: ',...
    num2str(user.MemAvailableAllArrays./1e6),'MB',newline,...
    'MemUsedMATLAB: ',num2str(user.MemUsedMATLAB./1e6),'MB']);
disp(' ');