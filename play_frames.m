function play_frames(mov,varargin)
% Author: Joshua Yeh
% Date created: 2017/11/16
%
%% SYNTAX
% play_frames(mov)
% play_frames(mov,pause_time)
% play_frames(moc,pause_time,flag)
% 
%% DESCRIPTION
% The purpose of this function is to accept a mov structure variable
% (output from the extract_frame function) and browse through the frames in
% a figure window. This function also allows the user to graphically browse
% through the frames by interacting with push buttons and edit boxes (very
% simple GUI).
% 
%% INPUT VARIABLES
% mov: a structure variable (output from the extract_frame function)
% pause_time: the amount of time (s) to pause between displaying each frame
% flag: a flag variables that determines whether to play through the frames
% when the function is running (1 = play, 0 = don't play).
%%
% Figure out the number of inputs
switch nargin
    case 1
        pause_time=0;
        flag=0;%don't play through frames
    case 2
        pause_time=varargin{1};
        flag=1;%play through frames
end

% Create a figure window
handles=guidata(figure(100));
handles.f1.f=gcf;clf(handles.f1.f);
handles.f1.s1=axes;
handles.virgin_flag=1;

%Show each frame
for dum=1:size(mov,2)    
    imagesc(handles.f1.s1,mov(dum).CData);
    axis image;
    
    %Include absolute frame number and timestamp
    title(handles.f1.s1,['Frame (abs): ',num2str(mov(dum).abs_frame_index),'/',num2str(mov(end).abs_frame_index(end)),...
        ', Timestamp: ',num2str(mov(dum).CurrentTime)]);
    drawnow;%draw figure graphics
    pause(pause_time);
    
    %Break out of for loop if flag is turned off
    if flag==0
        break
    end
end

% Create uicontrol
handles.prev=uicontrol('string','<','units','normalized',...
    'tooltipstring','Previous frame','position',[0.05 0.01 0.05 0.05],...
    'callback',@prev_callback,'style','pushbutton');
handles.next=uicontrol('string','>','units','normalized',...
    'tooltipstring','Next frame','position',[0.11 0.01 0.05 0.05],...
    'callback',@next_callback,'style','pushbutton');
handles.go=uicontrol('style','edit','units','normalized',...
    'tooltipstring','Go to frame','position',[0.17 0.01 0.075 0.05],...
    'callback',@go_callback,'string',mov(dum).abs_frame_index);
handles.mov=mov;%store mov var structure into the handles structure
handles.current=mov(dum).abs_frame_index;%store the current frame number
handles.xlim0=handles.f1.s1.XLim;
handles.ylim0=handles.f1.s1.YLim;
guidata(handles.f1.f,handles);%create handles structure for the figure window

% This is the callback fcn for handles.prev pushbutton
    function prev_callback(hObject,eventdata)
        f=ancestor(hObject,'figure');
        handles=guidata(f);
        if handles.current==handles.mov(1).abs_frame_index
            disp('Cannot go back!');
            return
        end
        
        % get the xlim and ylim values
        xlim0=handles.f1.s1.XLim;
        ylim0=handles.f1.s1.YLim;
        
        cla(handles.f1.s1);
        abs_frame_index=[handles.mov(:).abs_frame_index];%abs frame index
        ii=find(abs_frame_index==handles.current);
        handles.current=abs_frame_index(ii-1);%go back one frame        
        
        %show frame
        imagesc(handles.f1.s1,handles.mov(ii-1).CData);        
        
        % set the xlim and ylim to original values
        if handles.virgin_flag==0
            axis image;
            set(handles.f1.s1,'xlim',xlim0,'ylim',ylim0);
        else
            axis image;
        end

        title(handles.f1.s1,['Frame (abs): ',num2str(abs_frame_index(ii-1)),'/',num2str(abs_frame_index(end)),...
        ', Timestamp: ',num2str(handles.mov(ii-1).CurrentTime)]);
        drawnow;
        handles.go.String=handles.current;
        handles.virgin_flag=0;
        guidata(handles.f1.f,handles);%update the handles structure
    

% This is the callback fcn for handles.next pushbutton
    function next_callback(hObject,eventdata)
        f=ancestor(hObject,'figure');
        handles=guidata(f);
        if handles.current==handles.mov(end).abs_frame_index
            disp('Cannot go forward!');
            return
        end
        
        % get the xlim and ylim values
        xlim0=handles.f1.s1.XLim;
        ylim0=handles.f1.s1.YLim;
        
        cla(handles.f1.s1);
        abs_frame_index=[handles.mov(:).abs_frame_index];%abs frame index
        ii=find(abs_frame_index==handles.current);
        handles.current=abs_frame_index(ii+1);%go forward one frame        

        %show frame
        imagesc(handles.f1.s1,handles.mov(ii+1).CData);        
        
        % set the xlim and ylim to original values
        if handles.virgin_flag==0
            axis image;
            set(handles.f1.s1,'xlim',xlim0,'ylim',ylim0);
        else
            axis image;
        end
        
        title(handles.f1.s1,['Frame (abs): ',num2str(abs_frame_index(ii+1)),'/',num2str(abs_frame_index(end)),...
        ', Timestamp: ',num2str(handles.mov(ii+1).CurrentTime)]);
        drawnow;
        handles.go.String=handles.current;
        handles.virgin_flag=0;
        guidata(handles.f1.f,handles);%update the handles structure
    

% This is the callback fcn for handles.go edit text box
    function go_callback(hObject,eventdata)
        f=ancestor(hObject,'figure');%get figure parent
        handles=guidata(f);%get handles structure from figure
        n=str2double(hObject.String);%get frame number
        if n>handles.mov(end).abs_frame_index||n<handles.mov(1).abs_frame_index
            disp('Frame # outside of frame range!');
            return
        end
        
        % get the xlim and ylim values
        xlim0=handles.f1.s1.XLim;
        ylim0=handles.f1.s1.YLim;
        
        cla(handles.f1.s1);
        abs_frame_index=[handles.mov(:).abs_frame_index];%abs frame index
        handles.current=n;%update current frame
        
        %show frame
        ii=find(abs_frame_index==handles.current);%locate frame index
        imagesc(handles.f1.s1,handles.mov(ii).CData);        
        
        % set the xlim and ylim to original values
        if handles.virgin_flag==0
            axis image;
            set(handles.f1.s1,'xlim',xlim0,'ylim',ylim0);
        else
            axis image;
        end
        
        title(handles.f1.s1,['Frame (abs): ',num2str(abs_frame_index(ii)),'/',num2str(abs_frame_index(end)),...
        ', Timestamp: ',num2str(handles.mov(ii).CurrentTime)]);
        drawnow;
        handles.virgin_flag=0;
        guidata(handles.f1.f,handles);%update the handles structure
    