function [closed_frame, opened_frame]=extract_open_close(filenames,varargin)
% Joshua Yeh
% Date created: 19/04/05
%% DESCRIPTION
% The purpose of this function is to extract the open/close configurations
% of samples undergoing cyclic fatigue loading. The selected frames for
% analysis will fo through a secondary filter, so that the selected frames
% are not blurred.
%
%% INPUT PARAMETERS
% filenames: The name of the avi file containing the relevant frames to be
% analyzed (NOT IMPLEMENTED). This variables can also be a structure
% varibles containing a series of images to be analyzed.
%    -If the variable is a series of images, the structure must be formatted in
%    the same way as the output from the built-in MATLAB fcn, dir. See dir
%    for more information (doc dir).
%   -If the variable is a series of images, the analysis will process in
%   the order in which the filenames are listsed in the structure.
% 
% varargin: '<fieldname>',<value>
% 'ROI': region of interest encompassing the clamps, which will be used to
% determine maximum amplitude
% 
% 'rotate': rotate frames/images in deg. counterclockwise
% 
% 'thresh': threshold for binarization
% 
%% OUTPUT PARAMETERS
% 
%% PARSE THE INPUT VARIABLES

narginchk(1,inf);%check number of inputs is correct
params=inputParser;
params.CaseSensitive=false;
params.addParameter('ROI',[1 500 1 10],@(x) isnumeric(x));
params.addParameter('rotate',[],@(x) isnumeric(x));
params.addParameter('thresh',[],@(x) isnumeric(x));
params.parse(varargin{:});

% Extract out values from parsed input
ROI=params.Results.ROI;
rot=params.Results.rotate;
thresh=params.Results.thresh;

% Prepare time lines
f1=my_fig(1,{[3 1 1],[3 1 2],[3 1 3]});
set(f1.s1,'ydir','reverse');
set(f1.s3,'ycolor','none','box','off','ylim',[0 5]);
f1.p1=plot(f1.s3,nan,nan,'ko','displayname','all');
f1.p2=plot(f1.s3,nan,nan,'ro','displayname','close');
f1.p3=plot(f1.s3,nan,nan,'bo','displayname','open');
f1.p4=plot(f1.s3,nan,nan,'rx','linewidth',2,'markersize',8);
f1.p5=plot(f1.s1,nan,nan,'k-','tag','disp');
f1.p6=plot(f1.s2,nan,nan,'k-','tag','slope');
ylabel(f1.s1,'rel. pos. (px)');
ylabel(f1.s2,'local slope');
xlabel(f1.s3,'frame index');
f1.L=legend([f1.p1 f1.p2 f1.p3],'orientation','horizontal','location','north');
f1.L.FontSize=10;
linkaxes([f1.s1 f1.s2 f1.s3],'x');

% Add a custom datacursormode to the timeline figure
dcm_obj=datacursormode(f1.f);
set(dcm_obj,'UpdateFcn',{@show_frame,filenames});

% show the frame associated with the current/selected frame index
f2=my_fig(2);
f2.I=imagesc(f2.s1,nan);
f2.I.Tag='open_close_frame';
f2.p1=plot(f2.s1,nan,nan,'r-');
axis(f2.s1,'image');
set(f2.s1,'fontsize',8,'ydir','reverse');
center_axes(f2.s1,'margins',30);

% show the mean intensity profile
f3=my_fig(3);
f3.p1=plot(f3.s1,nan,nan,'k-');
f3.p2=plot(f3.s1,nan,nan,'r-');
xylabels(f3.s1,'array index','intensity (a.u.)');
center_axes(f3.s1);

% check to see if filenames is an avi or a series of images
if isstruct(filenames)
    format='images';
    N=1:size(filenames,1);%frame index array
    set(f1.p1,'xdata',N,'ydata',ones(size(N)));
    set(f1.s1,'xlim',[1 numel(N)]);
    disp('image series detected');
    disp([num2str(numel(N)),' frames detected']);
elseif ~isstruct(filenames)&&strcmp(filenames(end-3:end),'.avi')
    format='avi';
    disp('.avi video file detected');
end


%% Process the images

if strcmp(format,'images')%if filenames is a series of images
    % preallocate arrays
    kk=nan(size(filenames));
    slope=kk;
    for dum=1:size(filenames,1)
        filename=[filenames(dum).folder,'/',...
            filenames(dum).name];%full filename
        
        % Show the current image being analyzed
        set(f1.p4,'xdata',dum,'ydata',1);
        
        % Import the image
        I0=imread(filename);
        
        if ~isempty(rot)%chekc if need to rotate image
            I0=imrotate(I0,rot);
        end
        
        % extract out the subimage
        I1=I0(ROI(1):ROI(2),ROI(3):ROI(4),1);
        I2=nanmean(I1,2);% take mean in col dir
        if isempty(thresh)%define threshold value if user didnt define it
            thresh=graythresh(I2);
            thresh=thresh*max(I2);
        end
        
        % find index that intersects the threshold value
        kk(dum)=find(I2>=thresh,1,'first');
        
        % determine local slope around the intersection point
        slope(dum)=(I2(kk(dum)+1)-I2(kk(dum)-1))/2;
        
        set(f1.p5,'xdata',N,'ydata',kk);
        set(f1.p6,'xdata',N,'ydata',slope);
        
        % show ROI mean profile
        set(f3.p1,'ydata',I2,'xdata',1:numel(I2));
        set(f3.p2,'ydata',[thresh thresh],'xdata',[1 numel(I2)]);
        
        % show the raw image
        title(f2.s1,['frame ',num2str(dum)]);
        set(f2.I,'cdata',I0);
        set(f2.p1,'xdata',ROI([3 4 4 3 3]),'ydata',ROI([1 1 2 2 1]));
        uistack(f2.p1,'top');
        drawnow

    end
elseif strcmp(format,'avi')%if filenames is an avi video
    
end

% find peaks in rel. displacement
[pks,locs1]=findpeaks(kk);
[valleys,locs2]=findpeaks(-(kk-nanmean(kk)));
valleys=-valleys+nanmean(kk);

% find max slope
[~,idx_close]=max(slope(locs1));
[~,idx_open]=max(slope(locs2));

% update plots
set(f1.p2,'xdata',locs1,'ydata',ones(size(locs1)).*2);
set(f1.p3,'xdata',locs2,'ydata',ones(size(locs2)).*3);

plot(f1.s1,locs1,pks,'rx','displayname','close');
plot(f1.s1,locs2,valleys,'bx','displayname','open');
plot(f1.s1,locs1(idx_close),pks(idx_close),'rp','markerfacecolor','r');
plot(f1.s1,locs2(idx_open),valleys(idx_open),'bp','markerfacecolor','b');
plot(f1.s2,locs1,slope(locs1),'rx','displayname','close');
plot(f1.s2,locs2,slope(locs2),'bx','displayname','open');
plot(f1.s2,locs1(idx_close),slope(locs1(idx_close)),'rp','markerfacecolor','r');
plot(f1.s2,locs2(idx_open),slope(locs2(idx_open)),'bp','markerfacecolor','b');
plot(f1.s3,locs1(idx_close),2,'rp','markerfacecolor','r',...
    'displayname','best closed');
plot(f1.s3,locs2(idx_open),3,'bp','markerfacecolor','b',...
    'displayname','best opened');
set(findall(f1.f,'marker','p'),'markersize',12);

% output optimized open and close frame number
closed_frame=locs1(idx_close);
opened_frame=locs2(idx_open);

disp(['opened: ',num2str(opened_frame),'     closed: ',num2str(closed_frame)]);
disp(['fini: ',datestr(clock)]);

%% Useful fcns
function txt=show_frame(~,event_obj,filenames)

    % Customizes text of data tips
    pos = get(event_obj,'Position');
    type=event_obj.Target.DisplayName;
    
    % get rel. disp. data
    d1=findall(0,'tag','disp');
    disp_txt=['rel. pos.: ',num2str(d1.YData(pos(1)))];
    
    % get slope data
    d2=findall(0,'tag','slope');
    slope_txt=['slope: ',num2str(d2.YData(pos(1)))];
    
    txt = {['frame index: ',num2str(pos(1))],...
        type,disp_txt,slope_txt};
          
    % get handle of the image frame
    h=findall(0,'tag','open_close_frame');
    filename=[filenames(pos(1)).folder,'/',...
        filenames(pos(1)).name];
    I0=imread(filename);
    h.CData=I0;
    title(h.Parent,['frame ',num2str(pos(1))]);