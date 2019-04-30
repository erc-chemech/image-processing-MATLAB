function results=DA_integrate(filename,ffmat,dmap,varargin)
% Joshua Yeh
% 190315
%% DESCRIPTION
% The purpose of this script is to perform DA integration on a series of
% images. Note that this script is written only for straight crack edges.
% 
%% Input variables
% filename: filenames of the image that will be analyzed
% 
% ffmat: flat-field correction mat file
% 
% dmap: distortion correction mat file
% 
% NAME PAIR ARGUMENTS: DA_integrate(...'<fieldname>',<value>)
% 
% 'calfile': calibration mat file containing 'fo', a fit object that
% converts intenstiy to concentration (rel. %)
% 
% 'w2': pixel width when performing 2d-binning converting distortion
% corrected pixel coordinates-intensities
% 
% 'T1': threshold value that defines the region where stats are done based
% on the flat-field correction array
%
% 'res_p': length of the perpendicular line profile rel. to the edge in
% pxs
% 
% 'n_cuts': number of line cuts along the crack edge in pxs
% 
% 'thresh': background intensity subtraction (post correction)
% 
% 'binarize_thresh': threshold to binarize the image for edge detection
%
% 'r_divide': this controls the behavior of where to perform the
% integration in the image (um)
% 
% 'edge_detect': 'manual' or 'auto' (default) mode for defining the sample
% crack edge
% 
%% parse input variables
narginchk(3,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('calfile',[],@(x) ischar(x));
params.addParameter('w2',2,@(x) isnumeric(x));
params.addParameter('T1',0.4,@(x) isnumeric(x));
params.addParameter('res_p',150,@(x) isnumeric(x));
params.addParameter('n_cuts',150,@(x) isnumeric(x));
params.addParameter('thresh',188,@(x) isnumeric(x));
params.addParameter('binarize_thresh',1400,@(x) isnumeric(x));
params.addParameter('r_divide',500,@(x) isnumeric(x));
params.addParameter('edge_detect','auto',@(x) ischar(x)&&...
    (strcmp(x,'auto')||strcmp(x,'manual')));
params.parse(varargin{:});

% Extract out values from parsed input
calfile=params.Results.calfile;
w2=params.Results.w2;
T1=params.Results.T1;
res_p=params.Results.res_p;
n_cuts=params.Results.n_cuts;
thresh=params.Results.thresh;
binarize_thresh=params.Results.binarize_thresh;
r_divide=params.Results.r_divide;
edge_detect=params.Results.edge_detect;

% import image
I=import_tiff_stack(filename,1,'skip',1,'silence',1);

% define additional variables based on parsed inputs
H=w2*res_p;% length of line profile in um
LE=w2*n_cuts;% target length of edge

% Load flat-field correction map
load(ffmat,'IH');

% Load distortion map
load(dmap,'dx','dy');

% check to see if user provided calibration file
if ~isempty(calfile)
    load(calfile,'fo');
end

%% Prepare dynamic plots
f1=my_fig(1);
f1.c=colorbar(f1.s1);
axis(f1.s1,'image');
colormap(f1.s1,'parula');
xylabels(f1.s1,'x (\mum)','y (\mum)');
center_axes(f1.s1,'margins',10);

f2=my_fig(2);
colormap(f2.s1,'bone');
xylabels(f2.s1,'distance from edge (\mum)','integrated intensity');
center_axes(f2.s1,'margins',10);

%% Perform image corrections

% I=import_tiff_stack(filename,1,'skip',1,'silence',1);% import image
plane_initial=(I.tiff_stack)./IH;%flat field correction
plane_initial(IH(:)<T1)=0;%remove overcorrect areas, defined by T1

% Convert px coord to cartesian coord
[x,y]=meshgrid(1:size(plane_initial,2),1:size(plane_initial,1));
x=x.*I.res;%convert to real position in um
y=y.*I.res;
intensity=plane_initial(:);%intensity values

% Perform distortion correction
x=x-dx;
y=y-dy;

% Turn 2d array into column array
[x,y]=prepareCurveData(x,y);

% Performing 2D binning of current image
[I1,A,B]=coord2image(x,y,intensity,w2,'mean');
I1=I1-thresh;%remove background threshold
I1=medfilt2(I1,[3 3]);%2d median filter
I1(I1<0)=nan;%set values below 0 as nan

%% fit a linear line to straight edge

% Show the processed image
cla(f1.s1);
scatter(f1.s1,A(:),B(:),'.','cdata',I1(:));
set(f1.s1,'clim',[0 4095],'ydir','reverse');

% Binarize the image to prepare for linear edge fitting
x2=A(I1>binarize_thresh);%x coord of (um) intensities above threshold
y2=B(I1>binarize_thresh);%y coord of (um) intensities above threshold

% fit a linear line using total least squares
line_fit=polyfit(x2,y2,1);%first use linear least sq. as an initial guess
mypoly=@(x,p) polyval(p,x);
[~,line_fit,~]=numerFminS(mypoly,2,x2,y2,'initial',line_fit);%tot. least sq.
r2_fit=polyval(line_fit,x2);%resulting fitted y-coord (um) based on x2
c3=linspace(min(x2),max(x2),1000);%x pos of fitted edge (um) (higher res)
r2_fit2=polyval(line_fit,c3)';%corresp. y pos of fitted edge (higher res)
slope=line_fit(1);%slope of the fitted line

if strcmp(edge_detect,'auto')
    % show the edge line
    plot3(f1.s1,x2,r2_fit,ones(size(x2)).*1e5,'r-','linewidth',2);
elseif strcmp(edge_detect,'manual')
    drawnow;
    kk=find(x2>0&x2<nanmax(A(:))&r2_fit>0&r2_fit<nanmax(A(:)));
    edge_h=imline(f1.s1,[x2(kk(1)) r2_fit(kk(1));x2(kk(end)) r2_fit(kk(end))]);
    setColor(edge_h,'r');
    figure(f1.f);
    
    % Create dialogue button for manual edge line confirmation
    qf=figure;
    set(qf,'position',[500 678 200 100],'toolbar','none');
    tb=uicontrol(qf,'position',[1 30 198 40],'string',...
        'Confirm edge line','value',0,'style','togglebutton');
    
    % wait for user to confirm edge selection
    waitfor(tb,'value',1);
    delete(qf);%delete confirmation window
    
    % get user define edge line position, slope, and y-intercept
    pos=getPosition(edge_h);
    slope=(pos(1,2)-pos(2,2))/(pos(1,1)-pos(2,1));
    intercept=pos(1,2)-slope*pos(1,1);
    line_fit=[slope intercept];
    
    % recalculate edge line based on user defined line
    r2_fit=polyval(line_fit,x2);%resulting fitted y-coord (um) based on
    c3=linspace(min(pos(:,1)),max(pos(:,2)),1000);%x pos of fitted edge (um) (higher res)
    r2_fit2=polyval(line_fit,c3)';%corresp. y pos of fitted edge (higher res)
    
end

set(f1.s1,'xlim',[0 850],'ylim',[0 850],'clim',[10 4095],...
    'ydir','reverse');

% determine range along the edge to perform analysis
kk3=knnsearch(r2_fit2,r_divide);
c3_end=c3(kk3);%end x pos (um)
c3_start=c3_end-(LE)/sqrt(1+slope^2);%starting x pos (um)

% edge length (microns) (recalc. to make sure the math is right)
edge_L=sqrt((line_fit(1)*abs(c3_end-c3_start))^2+...
   (abs(c3_end-c3_start))^2);
if abs(edge_L-LE)>1e-6
    disp('something went wrong with the edge length calculations');
    disp(edge_L); %edge_L should be the same as LE
    disp(LE);
end

% cycle though along edge and get the intensity profile perpendicular to
% the fitted crack edge
improfiles=nan(n_cuts,res_p);%preallocate intensity profile line
count=1;%cycle number
cqs=linspace(c3_start,c3_end,n_cuts);

for cq=cqs
    % determine the perpendicular line
    perp_line_fit=[-1/line_fit(1),...
        polyval(line_fit,cq)+1/line_fit(1)*cq];

    phi=atand(perp_line_fit(1));%slope in degrees
    cq_end=H.*cosd(phi)+cq;%end x pos (um) of perpendicular line profile

    % extract and store improfile (intensity profile along line)
    xi=[cq cq_end];
    yi=polyval(perp_line_fit,xi);
    improfiles(count,:)=improfile(A,B,I1,xi,yi,res_p,'bicubic');
    
    % plot perpendicular line
    plot3(f1.s1,cq:cq_end,polyval(perp_line_fit,cq:cq_end),...
        ones(size(cq:cq_end)).*1e5,...
        'r-','linewidth',2);

    drawnow;
    count=count+1;
end

%% integrate image intensities

% Need to find the cutoff for integration
% This will be based on the cutoff noise, if the average intensity
% at a distance fro m the edge is less than then noise, then set the
% intensity as nan
improfiles2=improfiles;
improfiles2(isnan(improfiles2))=0;%convert nans to 0s
improfiles_mean=nanmean(improfiles2,1);%mean of all line intensity profiles
improfiles_std=nanstd(improfiles2,1,1);%std of all line intensity profiles
cut_i=find(improfiles_mean<thresh);% find cutoff to thresh noise
if ~isempty(cut_i)
    [~,idx]=nanmax(improfiles_mean);
    cut_i2=find(cut_i>idx,1,'first');
    cut_i=cut_i(cut_i2);
else%if all intensities are above thresh noise
    cut_i=numel(improfiles_mean);
end

Hs=linspace(0,H,size(improfiles2,2));%integration dist. from edge (um)
SI=cumsum(improfiles_mean(1:cut_i),2);%integrated intensity
if isempty(SI)
    disp('Something went wrong with the integration!');
end

%% store the results
results.Hs=Hs(1:cut_i);%integration dist. from edge (um)
results.SI=SI;%surface integrated intensity

% Statistical calc.
results.improfiles_avg=improfiles_mean(1:cut_i);%mean
results.improfiles_std=improfiles_std(1:cut_i);%stdev
results.improfiles=improfiles2;%intensity profiles for each perpend. line
results.params=params;%store the parameters of the analysis in the results

% check to see if user provided calibration file
if ~isempty(calfile)
    % Convert intensities to concentrations
    results.conc_avg=fo(results.improfiles_avg);
    results.conc_std=fo(results.improfiles_std);
end

% plot suf. integ. results
plot(f2.s1,Hs(1:cut_i),SI,'o-');
center_axes(f2.s1);

% plot the histogram intensities
f3=my_fig(3);
histogram(f3.s1,I.tiff_stack(:));
xylabels(f3.s1,'intensity (a.u.)','counts');
set(f3.s1,'yscale','log');
center_axes(f3.s1);

% plot avg intensity profile
f4=my_fig(4);
plot(f4.s1,Hs,improfiles_mean);
xylabels(f4.s1,'distance from crack edge (\mum)','intensity (a.u.)');
center_axes(f4.s1);


