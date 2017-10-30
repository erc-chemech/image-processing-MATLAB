function mov=extract_frames(filename,frame_range)
% This functions serves to extract a range of frames from a video file. The
% frames are stored in a structure variable.
% Author: Joshua Yeh
% Date created: 2017/10/28
%
% INPUT VARIABLES
% filename: name of video file (RGB file)
% frame_range: range of frames [first last] to extract from video file
%   -If user inputs the string 'all' for the frame_range value, the
%   function will import all of the frames from the video file.
% 
% OUTPUT VARIABLES
% mov: a structure variable containing the extracted frames
%   fieldnames:
%       CData: frame data
%       cmap: colormap
%       abs_frame_index: absolute frame index (in reference to video file)

% Turn off hardware acceleration to prevent VideoReader from crashing (an
% apparant graphics card issue that crashes MATLAB)

% Check to see if user inputs 'all' for abs_frame_index
if strcmp(frame_range,'all')==1&&ischar(frame_range)==1
    disp('This fcn will attempt to import all frames.');    
elseif strcmp(frame_range,'all')==0&&ischar(frame_range)==1
    disp('Did not recognize char input for abs_frame_index.');
    disp('Frame import aborted.');
    return
else
end

matlab.video.read.UseHardwareAcceleration('off')

%Create a VideoReader object
Vidobj = VideoReader(filename);

% Turn back on the hardware acceleration
matlab.video.read.UseHardwareAcceleration('on')

% Import frames from video into mov struct variable
mov = struct('CData',zeros(Vidobj.Height,Vidobj.Width,3,2,'uint8'),'cmap',[],...
    'abs_frame_index',[]);
n=1;%counter absolute frame index for the video file
k=1;%counter index for mov structure

disp('Importing frames (this can take awhile)');
while hasFrame(Vidobj)
    % If the frame falls within the designated frame range, store the frame
    % in the mov structure variable.
    if ischar(frame_range)==0%if user specifies range of frames to import
        if n>=frame_range(1)&&n<=frame_range(2)
            mov(k).CData=readFrame(Vidobj);%store the frame
            mov(k).abs_frame_index=n;%store the absolute frame index
            if mod(k,100)==0%update user every 100 frames are read
                disp([num2str(k),' frames imported']);
            end
            k=k+1;
        else % Otherwise, do not store the frame.
            readFrame(Vidobj);%need this to advance to the next frame
        end
    % If user specifies to import all frames.
    elseif ischar(frame_range)==1&&strcmp(frame_range,'all')
        mov(k).CData=readFrame(Vidobj);%store the frame
        mov(k).abs_frame_index=n;%store the absolute frame index
        if mod(k,100)==0%update user every 100 frames are read
            disp([num2str(k),' frames imported']);
        end
        k=k+1;
    end
    n=n+1;
end
disp('Import success');
disp(['Imported ',num2str(k-1),' frames']);