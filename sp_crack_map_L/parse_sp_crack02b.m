function parse_sp_crack02b(filename,frame_n,varargin)
% Author: Joshua Yeh
% Date created: 2018-04-18
% Updated on: 2018-05-28
% 
%% DESCRIPTION (FOR LOADING ANALYSIS)
% This script accepts a video frame and analyzes the pixels in the frame.
% The pixels are categorized as loading, background, or other.
% The pixels are plotted on a RCC vs. BCC chromatic plot.
% 
%% INPUT VARIABLES
% filename: the name of the file that will be used to perform the analysis
% 
% frame_n: the frame in which the analysis will be performed on
% 
% varargin: <'field_name'>,<'value'>
    %
    % 'skip': A toggle flag where a value of 1 signifies to the function to
    % skip manual import of the frame from the video file. A value of 0
    % signifies the opposite.
    %
    % 'ref_frame': the reference frame representing the beginning of the
    % experiment (default is 1)
    %
    % 'ROI': the sample region of interest in the frame being analyzed
    % 
    % 'ref_ROI': the sample region of interest in the reference frame
    %
    % 'white_ref_ROI': the region of interest representing the white area
    % in the reference frame
    %
    % 'background_thresh': the background threshold value that
    % distinguishes the background from the foreground (for masking
    % purposes) The threshold value is taken rel. to the total px
    % intensity.
    %
    % 'dark_thresh': the dark thereshold value that excludes the dark
    % regions in the sample (for masking purposes) The threshold value is
    % taken rel. to the total px intensity.
    %
    % 'stress_calc': toggle (1 or 0, 1 is default) for performing stress
    % association of identified loading and unloading responses
    %
    % 'export_var': export variables from the function workspace to the
    % base workspace
    %
    % 'RCC_thresh': percentile (0 to 1) threshold distinguishing signal
    % from noise
    %
    % 'BCC_thresh': percentile (0 to 1) threshold distinguishing sigal from
    % noise
    %
    % 'res': this parameter controls the hexagon size for the honeycomb
    % 2d histrogram plots
    % 
    % 'import_mask': this is a flag (0 (default) or 1) to manually import
    % an image mask and use it in addition to the masking algorithm written
    % in this function
    %
    % 'import_mask_file': this is the filename (be sure to include the
    % path) of the image mask the user can import and use
    %
    % 'noise_filter': this is a flag (0 (default) or 1) for applying a
    % noise filter after the picel parsing process
    %
    % 'video_mode': this is a flag (0 (default) or 1) for running the
    % algorithm in video mode
    %
    % 'initialize': determine whether ('normal', default) or not ('only')
    % to run the intial frame extraction sequence
    %
    % 'LoadFitFile': the .mat file containing the fit objects of the
    % loading curve
    % 
    % 'CSMFile':the .fig containing the color-stress map
    %
    % 'res': the step size for determining the x and y edges for the
    % hexagonal histogram binning and plotting for the chromatic density
    % plots
    % 
    % 'max': the maximum counts at saturation for the honeycomb 2d
    % histogram
    %
    % 'Ires': image resolution in px2mm
    %
    % 'map_value': stress in value (scalar) ('stress') or the energy
    % density ('energy_density')
    %
    % 'calc_true_stress': toggle for calc. the true stress map (do not
    % calc. 0    or    calc true stress 1)
    %
    % 'true_stress_mat': mat file associated with calculating the true
    % stress. This mat file must contain a variable named 'TS_NS', which is
    % a fit object correlating true stress vs. nominal stress.
    %
    % 'corr_factor': a correction factor that takes into account the diff.
    % distance between the camera between the calibration and fracture exp/
    %
    % 'S2E_conversion': file that converts stresses to energy density
    %
    % 'm_file': mechanical data (xls file) generated recalc_lambda fcn
    %
    % 'rotate': rotate frames in degrees
    % 
%% Parse input variables

narginchk(1,inf);%check the number of input values
params=inputParser;
params.CaseSensitive=false;
params.addParameter('skip',1,@(x) islogical(x)|x==1|x==0);
params.addParameter('ref_frame',1,@(x) isnumeric(x));
params.addParameter('ROI',[500 800 60 270],@(x) isnumeric(x)&numel(x)==4);
params.addParameter('ref_ROI',[800 1000 250 300],...
    @(x) isnumeric(x)&numel(x)==4);
params.addParameter('white_ref_ROI',[800 1000 500 600],...
    @(x) isnumeric(x)&numel(x)==4);
params.addParameter('background_thresh',565,@(x) isnumeric(x));
params.addParameter('dark_thresh',340,@(x) isnumeric(x));
params.addParameter('stress_calc',1,@(x) islogical(x)|x==1|x==0);
params.addParameter('export_var',1,@(x) islogical(x)|x==1|x==0);
params.addParameter('RCC_thresh',0.05,@(x) isnumeric(x)&x>0&x<1);
params.addParameter('BCC_thresh',0.95,@(x) isnumeric(x)&x>0&x<1);
params.addParameter('import_mask',0,@(x) isnumeric(x));
params.addParameter('import_mask_file',' ',@(x) ischar(x));
params.addParameter('noise_filter',0,@(x) isnumeric(x));
params.addParameter('video_mode',0,@(x) isnumeric(x));
params.addParameter('initialize','normal',...
    @(x) ischar(x)&(strcmp(x,'only')|strcmp(x,'normal')));
params.addParameter('LoadFitFile',...
    '../cyclic loading results/TN135_loading_fits.mat',...
    @(x) ischar(x));
params.addParameter('CSMFile','CSM_TN135.fig', @(x) ischar(x));
params.addParameter('res',1e-3,@(x) isnumeric(x));
params.addParameter('max',100,@(x) isnumeric(x)&x>0);
params.addParameter('Ires',38.70,@(x) isnumeric(x));
params.addParameter('S2E_conversion','../G matching/TN135_E.mat');
params.addParameter('map_value','stress');
params.addParameter('calc_true_stress',0,@(x) islogical(x)|x==1|x==0);
params.addParameter('true_stress_mat','',@(x) ischar(x));
params.addParameter('corr_factor',1,@(x) isnumeric(x));
params.addParameter('m_file',[],@(x) ischar(x));
params.addParameter('rotate',0,@(x) isnumeric(x));
params.parse(varargin{:});

%Extract out values from parsed input
skip=params.Results.skip;
ref_frame=params.Results.ref_frame;
ROI=params.Results.ROI;
ref_ROI=params.Results.ref_ROI;
white_ref_ROI=params.Results.white_ref_ROI;
background_thresh=params.Results.background_thresh;
dark_thresh=params.Results.dark_thresh;
stress_calc=params.Results.stress_calc;
export_var=params.Results.export_var;
import_mask=params.Results.import_mask;
noise_filter=params.Results.noise_filter;
video_mode=params.Results.video_mode;
init=params.Results.initialize;
LoadFitFile=params.Results.LoadFitFile;
CSMFile=params.Results.CSMFile;
BCC_thresh=params.Results.BCC_thresh;
RCC_thresh=params.Results.RCC_thresh;
res=params.Results.res;
max1=params.Results.max;
Ires=params.Results.Ires;
S2E_conversion=params.Results.S2E_conversion;
map_value=params.Results.map_value;
calc_ts=params.Results.calc_true_stress;
ts_mat=params.Results.true_stress_mat;
corr_factor=params.Results.corr_factor;
m_file=params.Results.m_file;
rot=params.Results.rotate;

if video_mode==1
    visible='off';
else
    visible='on';
end

% extract stress strain data
if ~isempty(m_file)
    [num,txt,raw]=xlsread(m_file);
    
    %get stress value
    k=find(strcmp(txt,'Stress (MPa)'));
    [kr,kc]=ind2sub(size(txt),k);
    m_stress=[raw{kr+1:end,kc}]';
    
    %get time value
    k=find(strcmp(txt,'Time (s)'));
    [kr,kc]=ind2sub(size(txt),k);
    m_time=[raw{kr+1:end,kc}]';
    
    %get lambda value
    k=find(strcmp(txt,'Lambda'));
    [kr,kc]=ind2sub(size(txt),k);
    m_lam=[raw{kr+1:end,kc}]';
    
    %get abs frame number
    k=find(strcmp(txt,'frame #'));
    [kr,kc]=ind2sub(size(txt),k);
    frame_idx=[raw{kr+1:end,kc}]';
    ii_fn=find(frame_idx==frame_n);
end


%% LOAD VIDEO AND FITS
disp(['fcn performing on frame: ',num2str(frame_n)]);
[~,name,~]=fileparts(filename);
if skip==0
    frame=extract_frames(filename,[ref_frame frame_n]);
    save(['L_',name,'_',num2str(frame_n),'.mat'],'frame');
else
    try
        load(['L_',name,'_',num2str(frame_n),'.mat']);
    catch
        disp('mat file was not found, attempting to extract from file');
        frame=extract_frames(filename,[ref_frame frame_n]);
        save(['L_',name,'_',num2str(frame_n),'.mat'],'frame');
    end
end

% Load the loading fit curves for TN135
LFF=load(LoadFitFile,'RCCvBCC','stressvBCC');

disp('Frame loaded!');

if strcmp(init,'only')
    disp(['Fini: ',datestr(clock)]);
    return
end

%% PERFORM COLOR CORRECTIONS

% Check to see if the frame needs to be rotated
if rot~=0
    frame01=imrotate(frame(1).CData,rot);
    frame02=imrotate(frame(2).CData,rot);
else
    frame01=frame(1).CData;
    frame02=frame(2).CData;
end

% Crop our ROI
subimage=frame02(ROI(1):ROI(2),ROI(3):ROI(4),:);
ref=frame01(ref_ROI(1):ref_ROI(2),ref_ROI(3):ref_ROI(4),:);

% Define white reference area
white_ref=frame01(white_ref_ROI(1):white_ref_ROI(2),...
    white_ref_ROI(3):white_ref_ROI(4),:);

% Apply color correction
frame_corr=rgb_correction(subimage,white_ref,'simple',200);%for analysis
frame_corr_sum=sum(frame_corr,3);
frame_corr2=rgb_correction(subimage,white_ref,'simple',200);%for display
ref_corr=rgb_correction(ref,white_ref,'simple',200);%for analysis
ref_corr_sum=sum(ref_corr,3);


%% PERFORM FRAME MASKING

% Convert frame_corr to intensity matrix
I1=sum(frame_corr,3);

% Remove the background
ii=find(I1>background_thresh);
template=zeros(size(frame_corr(:,:,1)));% create a template of zeros
template(ii)=1;%based on the background threshold, set background as 1s
template=bwareaopen(template,500);%remove small isolated regions of 1s
template=bwmorph(template,'thicken',10);%thicken the main background area
template=bwareaopen(-template+1,5000);%remove small isolated regions of 0s
ii1=find(template==0);
frame_corr=clean(frame_corr,ii1);

% Remove dark regions attributed to dust and edges
ii=find(I1<dark_thresh);
frame_corr=clean(frame_corr,ii);

template=zeros(size(frame_corr(:,:,1)));%create template of 0s
template(find(~isnan(frame_corr(:,:,1))))=1;%set elements that are #s as 1s
template=bwareaopen(template,5000);%remove isolated areas of 1s
ii2=find(template==0);
frame_corr=clean(frame_corr,ii2);
    
% apply any additional imported masks specified by the user
if import_mask==1
    import_mask_file=params.Results.import_mask_file;
    imported_mask=imread(import_mask_file);
    if size(imported_mask,3)>1
        imported_mask=rgb2gray(imported_mask);
    end
    imported_mask=imbinarize(imported_mask);%binarize image
    template=template.*imported_mask;
    ii2=find(template==0);
    frame_corr=clean(frame_corr,ii2);
end

%% PERFORM COLOR CALCULATIONS

% Calculate RGB ratios
[R_ratio,G_ratio,B_ratio]=rgb_ratio(frame_corr);
[R_ratio0,G_ratio0,B_ratio0]=rgb_ratio(ref_corr);

% Calculate chromatic change
% main frame
RCC=(R_ratio-nanmean(R_ratio0(:))).*corr_factor;
GCC=(G_ratio-nanmean(G_ratio0(:))).*corr_factor;
BCC=(B_ratio-nanmean(B_ratio0(:))).*corr_factor;

% reference frame
RCC0=R_ratio0-nanmean(R_ratio0(:));
GCC0=G_ratio0-nanmean(G_ratio0(:));
BCC0=B_ratio0-nanmean(B_ratio0(:));


%% PARSE PIXELS AS LOADING, OTHER, OR BACKGROUND

% import the 0-iso stress line from color-stress map for TN135
try delete(fcsm.f); catch; end
fcsm.f=openfig(CSMFile,'invisible');
fcsm.xdata=get(findall(fcsm.f,'type','surface'),'xdata');%BCC
fcsm.ydata=get(findall(fcsm.f,'type','surface'),'ydata');%RCC
try delete(fcsm.f); catch; end

% Create indexed array of BCC and RCC
disp('Parsing color coordinates');
% IA col label: px vector index,BCC,RCC,fitted RCC,
% flag for tagging pxs as loading etc., px row, px col
IA=[[1:numel(RCC(:))]',BCC(:),RCC(:),polyval(LFF.RCCvBCC,BCC(:)),...
    zeros(numel(BCC(:)),1),...%flag col
    ones(numel(BCC(:)),1),ones(numel(BCC(:)),1)];
thresh_BCC=quantile(BCC0(:),BCC_thresh);
thresh_RCC=quantile(RCC0(:),RCC_thresh);

% Determine the closest distance between the current pixel chromatic
% coordinate and the chromatic loading curve.
BCC_load=linspace(0.001,0.1,1e3)';
RCC_load=polyval(LFF.RCCvBCC,BCC_load);
[d_index,d]=knnsearch([BCC_load RCC_load],IA(:,2:3));

for dum=1:numel(RCC)
    [r,c]=ind2sub(size(BCC),dum);
    IA(dum,6:7)=[r,c];
    
    if isnan(IA(dum,2))||isnan(IA(dum,3))
        continue
    end
    
    % conditions for being part of the sample background
    if IA(dum,2)<=thresh_BCC && IA(dum,3)>=thresh_RCC
        IA(dum,5)=2;
    
    % conditions for being part of the loading behavior
    elseif d(dum)<sqrt(thresh_BCC^2+thresh_RCC^2)
        IA(dum,5)=3;
    end
    
    if mod(dum,1e4)==0
        disp([num2str(dum),' of ',num2str(numel(BCC)),' completed']);
    end
end

% Extract the indices associated with each category
m3=find(IA(:,5)==3);%indices associated with loading color change
m2=find(IA(:,5)==2);%indices associated with sample background
m0=find(IA(:,5)==0);%indices not associated with noise/loading color change
nan1=find(isnan(IA(:,5)));%indices associated with background

% NOISE FILTERING
if noise_filter==1%if flag is turned on, apply a noise filter to the parsed pxs

    % Apply noise filter on a smaller local scale (loading) This removes
    % isolated pixels.
    template=zeros(size(BCC));
    template(m3)=1;%set elements corresponding to unloading to 1
    local_sum=conv2(template,ones(3),'same');%compute local sum
    ii=find(local_sum<=2&local_sum>0);
    IA(ii,5)=0;
    % Refresh loading points
    m3=find(IA(:,5)==3);%indices associated with loading color change

    % Apply noise filter on a larger local scale
    template=zeros(size(BCC));
    template(m3)=1;%set elements corresponding to unloading to 1
    local_sum=conv2(template,ones(8),'same');%compute local sum
    ii=find(local_sum<=5&local_sum>0);
    IA(ii,5)=0;
    % Refresh loading points
    m3=find(IA(:,5)==3);%indices associated with loading color change
end

%% ASSOCIATE UNLOADING AND LOADING POINTS WITH STRESS VALUES

if stress_calc==1

    disp('Associating loading points to stresses');
    m3_s=nan(numel(m3),1);%var. containing nominal stress assoc. w/ loading
    m3_E=nan(numel(m3),1);%var. containing energy density values
    
    if calc_ts==1%if toggle for calc. true stress is turned on
        try
            load(ts_mat,'TS_NS');
            m3_ts=nan(numel(m3),1);%var. containing true stress assos. w/ load.
        catch
            calc_ts=0;
            disp(['ts_mat file is not found or does not'...
                ' contain the ''TS_NS'' fit object']);
            disp('calc_ts toggle variable is set to 0');
        end
    end
    
    % load mat file needed to convert nominal stress to energy density
    S2E=load(S2E_conversion);
    
    % process loading points
    for dum=1:numel(m3)
        ii=d_index(m3(dum));%index of BCC pt closest to loading curve
        m3_s(dum)=polyval(LFF.stressvBCC,BCC_load(ii));%nominal stress in MPa
        m3_E(dum)=spline(S2E.stress,S2E.E,m3_s(dum));%energy density in MJ/m^3
        
        if calc_ts==1%det. true stress (MPa) based on correl. w/ nominal stress
            m3_ts(dum)=TS_NS(m3_s(dum));%true stress in MPa
        end

        if mod(dum,1e4)==0%update user on status
            disp([num2str(dum),' of ',num2str(numel(m3))]);
        end
    end
    disp('Finished associating loading points to stress');
else
    disp('Skipped stress association for loading and unloading points');
end

disp('Processing finished');

%% PLOT THE RESULTS

% Declare figure handles




% Figure 3


% Figure 4


% Figure 5


% Figure 6


% Figure 7


% Figure 8

% xylabels(f8.s1,'Blue chromatic change','Red chromatic change');

%% %%%%%%%%%%%%%%%%%%% figure 1 %%%%%%%%%%%%%%%%%%%%%

f1=my_fig(1,{[1 2 1] [1 2 2]},'visible',visible);
axis([f1.s1 f1.s2],'image');
set(f1.f,'name','Raw extracted frame');

imagesc(f1.s1,uint8(ref_corr));
imagesc(f1.s2,uint8(frame_corr));

%% %%%%%%%%%%%%%%%%%%% figure 2 %%%%%%%%%%%%%%%%%%%%%

f2=my_fig(2,{[3 2 1] [3 2 2] [3 2 3] [3 2 4] [3 2 5] [3 2 6]},...
    'fontsize',10,'visible',visible);
set(f2.f,'name','Corrected analyzed frame');
title(f2.s1,'ref_{corr}');
title(f2.s2,'frame_{corr}');
title(f2.s3,'BCC and RCC for ref');
title(f2.s4,'BCC and RCC for frame');
title(f2.s5,'ref intensities');
title(f2.s6,'frame intensities');
set(f2.s3,'xlim',[-0.03 0.03]);
set(f2.s4,'xlim',[-0.03 0.03]);
edges=0:2:256;
histogram(f2.s1,ref_corr(:,:,1),edges,'facecolor','r','facealpha',0.5);
histogram(f2.s1,ref_corr(:,:,2),edges,'facecolor','g','facealpha',0.5);
histogram(f2.s1,ref_corr(:,:,3),edges,'facecolor','b','facealpha',0.5);
histogram(f2.s2,frame_corr2(:,:,1),edges,'facecolor','r','facealpha',0.5);
histogram(f2.s2,frame_corr2(:,:,2),edges,'facecolor','g','facealpha',0.5);
histogram(f2.s2,frame_corr2(:,:,3),edges,'facecolor','b','facealpha',0.5);
histogram(f2.s3,BCC0(:),'facecolor','b','facealpha',0.5);
histogram(f2.s3,RCC0(:),'facecolor','r','facealpha',0.5);
histogram(f2.s3,GCC0(:),'facecolor','g','facealpha',0.5);
histogram(f2.s4,BCC(:),'facecolor','b','facealpha',0.5);
histogram(f2.s4,RCC(:),'facecolor','r','facealpha',0.5);
histogram(f2.s4,GCC(:),'facecolor','g','facealpha',0.5);
histogram(f2.s5,ref_corr_sum(:));
histogram(f2.s6,frame_corr_sum(:));

linkaxes([f2.s1 f2.s2],'x');
linkaxes([f2.s3 f2.s4],'x');
linkaxes([f2.s5 f2.s6],'x');

%% %%%%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%

f3=my_fig(3,{[1 1 1]},'visible',visible);
f3.f.Position(4)=565;
f3.s2=axes;
set(f3.s2,'position',[0.8 0.7 0.1 0.3]);
set(f3.f,'name','Color mapping of reference sample area');
xylabels(f3.s1,'blue chromatic change','red chromatic change');

plot3(f3.s1,BCC_load,RCC_load,...
    ones(1,numel(BCC_load)).*1000,'k--');
plot3(f3.s1,...
    [0 thresh_BCC thresh_BCC],...
    [thresh_RCC thresh_RCC 0],...
    ones(1,3).*1000,'r-','linewidth',2);

%transform chromatic datapoints to 2d histogram array
f3.surf1=my_honeycomb(f3.s1,BCC0(:),RCC0(:),'res',res,...
    'dn','reference','max',max1,'xEdges',[-0.1 0.1],...
    'yEdges',[-0.1 0.1]);
center_axes(f3.s1);
f3.c3=alphaColorbar(f3.s1,f3.surf1,'clab','Counts');

% image frame of reference area
imagesc(f3.s2,uint8(ref_corr));
axis(f3.s2,'image');
set(f3.s2,'box','off','xcolor','none','ycolor','none','xtick',[],...
    'ytick',[]);

axis(f3.s1,'image');
set(f3.s1,'xlim',[0 0.01],'ylim',[-0.01 0]);

%% %%%%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%

f4=my_fig(4,{[1 1 1]},'visible',visible);
f4.f.Position(4)=565;
f4.f.Alphamap=linspace(0,1,256);
set(f4.f,'name','Color mapping of activated sample area');
xylabels(f4.s1,'Blue chromatic change','Red chromatic change');
set(f4.s1,'xlim',[0 0.1],'ylim',[-0.1 0]);

plot3(f4.s1,BCC_load,RCC_load,...
    ones(1,numel(BCC_load)).*1000,'k--');

% Transform chromatic datapoints to 2d histogram array (m2)
f4.N_m3=my_honeycomb(f4.s1,BCC(m2),RCC(m2),'max',max1,'res',res,...
    'dn','background');

% % Transform chromatic datapoints to 2d histogram array (m0)
f4.N_m3=my_honeycomb(f4.s1,BCC(m0),RCC(m0),'max',max1,'res',res,...
    'dn','background');

% Transform chromatic datapoints to 2d histogram array (m3)
f4.N_m3=my_honeycomb(f4.s1,BCC(m3),RCC(m3),'max',max1,...
    'color',[255 186 0]./255,'res',res,'dn','loading');

f4.s1.Position=[0.2939 0.4407 0.5561 0.5327];
f4.c=alphaColorbar_stacked(f4.s1,max1,{'k',[255 186 0]./255});
f4.c.Position([2 4])=[0.13 0.1];
f4.c.YTick=[0.5 1.5];
f4.c.YTickLabel={'Background','Loading'};
f4.c.Box='off';
f4.c.TickDir='out';

%% %%%%%%%%%%%%%%%%%%% figure 5 %%%%%%%%%%%%%%%%%%%%%

f5=my_fig(5,{[1 1 1]},'marg_w',[0.01 0.01],'marg_h',[0.01 0.01],...
    'visible',visible);
f5.f.Alphamap=linspace(0,1,256);
f5.f.Position(3:4)=[667 645];

imagesc(f5.s1,uint8(rgb_correction(subimage,white_ref,'simple',255)));
plot(f5.s1,IA(m3,7),IA(m3,6),'b.');
patch(f5.s1,[0 78.94 78.94 0 0],[5 5 25 25 5],'w','edgecolor','none',...
    'facealpha',0.8);
plot(f5.s1,[0 48.94]+15,[10 10],'k-','linewidth',4);
text(f5.s1,32.25+10,10,'1 mm','horizontalalignment','center',...
    'verticalalignment','bottom');
axis(f5.s1,'image');
set(f5.s1,'box','off','xtick',[],'ytick',[],'xcolor','none',...
    'ycolor','none');

%% %%%%%%%%%%%%%%%%%%% figure 6 %%%%%%%%%%%%%%%%%%%%%

f6=my_fig(6);
f6.f.Position([3,4])=[340 330];
f6.s1.Layer='top';

if stress_calc==1
    copyobj(findall(f5.s1,'type','image'),f6.s1);% copy image frame
    axis(f6.s1,'image');
    f6.s1.Position(1)=0.07;
    set(f6.s1,'xtick',[],'ytick',[])

    % Plot chromatic coordinates on the image frame
    if strcmp(map_value,'stress')
        if calc_ts==0
            m3_m=m3_s;
            climm=[2 6];
            strm='nominal stress (MPa)';
        elseif calc_ts==1
            m3_m=m3_ts;
            climm=[2 8];
            strm='true stress (MPa)';
        end
    elseif strcmp(map_value,'energy_density')
        m3_m=m3_E;
        climm=[0.5 1.5];
        strm='energy density (MPa)';
    end
    
    zdata=px_coord2image(IA(m3,6),IA(m3,7),m3_m,size(ref_corr));
    far_field=nanmean(nanmean(zdata(:,end-10:end)));
    
    if ~isempty(m_file)
        notched_stress=m_stress(ii_fn);
        disp(['far field stress (MPa): ',num2str(far_field)]);
        disp(['notched stress (MPa): ',num2str(notched_stress)]);
        try
            ssfit=load(LoadFitFile,'ssfit');
            disp(['unnotched stress (MPa): ',num2str(ssfit.ssfit(m_lam(ii_fn)))]);
        catch
            disp('unable to load unnotched stress value');
        end
    end
    surf(f6.s1,zdata);
    cmap=jet(256);
    colormap(f6.s1,cmap);
    f6.s1.Position([2 4])=[0.4 0.55];
    set(f6.s1,'clim',climm);
    
    % Plot scalebar
    plot(f6.s1,[10 10+Ires],[10 10],'k-','linewidth',2);
    text(f6.s1,10+Ires/2,12,'1 mm','verticalalignment','bottom',...
        'horizontalalignment','center','fontname','Helvetica',...
        'fontsize',14);
    
    if ~isempty(m_file)
        text(f6.s1,10+Ires/2,24,['\lambda = ',num2str(m_lam(ii_fn),2)],...
            'verticalalignment','bottom',...
            'horizontalalignment','center','fontname','Helvetica',...
            'fontsize',14);
    %     text(f6.s1,10+Ires/2,36,['\sigma_n = ',num2str(m_stress(ii_fn),2),' MPa'],...
    %         'verticalalignment','bottom',...
    %         'horizontalalignment','center','fontname','Helvetica',...
    %         'fontsize',14);
    end
    
    center_axes(f6.s1,'margins',5);
end

%% %%%%%%%%%%%%%%%%%%% figure 7 %%%%%%%%%%%%%%%%%%%%%

if stress_calc==1
    
    f7=my_fig(7,{[1 2 1] [1 2 2]});
    f7.f.Alphamap=linspace(0,1,256);
    f7.f.Position(3:4)=[1100 660];
    f7.s1.Layer='top';
    xylabels(f7.s2,'Blue chromatic change','Red chromatic change');
    
    copyobj(findall(f5.s1,'type','image'),f7.s1);% copy image frame
    axis(f7.s1,'image');
    f7.s1.Position(1)=0.07;
    set(f7.s1,'xtick',[],'ytick',[])
    
    % Plot chromatic coordinates on the image frame
    scatter3(f7.s1,IA(m3,7),IA(m3,6),m3_m,'filled','cdata',m3_m,'sizedata',5);
    cmap=jet(256);
    colormap(f7.s1,cmap);
    set(f7.s1,'clim',climm);
    f7.s1.Position([2 4])=[0.4 0.55];
    
    % loading colorbar
    f7.c1=colorbar(f7.s1,'horizontal');
    
    f7.c1.Position(1:4)=...
        [f7.s1.InnerPosition(1) 0.291,...
        f7.s1.InnerPosition(3) 0.05];
    f7.c1.YLabel.String=strm;
    f7.c1.YLabel.FontSize=f7.s2.FontSize;
    set(f7.c1,'box','off','tickdir','both','limits',climm);
    
    % Plot scalebar
    plot(f7.s1,[10 10+Ires],[10 10],'k-','linewidth',2);
    text(f7.s1,10+Ires/2,12,'1 mm','verticalalignment','bottom',...
        'horizontalalignment','center','fontname','Helvetica',...
        'fontsize',14);
    
    if ~isempty(m_file)
        text(f7.s1,10+Ires/2,24,['\lambda = ',num2str(m_lam(ii_fn),2)],...
            'verticalalignment','bottom',...
            'horizontalalignment','center','fontname','Helvetica',...
            'fontsize',14);
    %     text(f6.s1,10+Ires/2,36,['\sigma_n = ',num2str(m_stress(ii_fn),2),' MPa'],...
    %         'verticalalignment','bottom',...
    %         'horizontalalignment','center','fontname','Helvetica',...
    %         'fontsize',14);
    end
    
    % chromatic-stress mapping
    copyobj(get(f4.s1,'children'),f7.s2);
    f7.s2.Position=[0.6047 0.4125 0.3447 0.5508];
    set(f7.s2,'xlim',[0 0.1],'ylim',[-0.1 0]);
    f7.c2=copyobj(f4.c,f7.f);
    f7.c2.Position=[0.605 0.125 0.345 0.1];
    f7.c2.YTick=[0.5 1.5];
    f7.c2.YTickLabel={'Background','Loading'};
    f7.c2.Box='off';
    f7.c2.TickDir='out';
    
    set(f7.s2,'xlim',[0 0.1],'ylim',[-0.1 0]);
    axis(f7.s2,'image');
    
end

% Export function workspace to base workspace
if export_var==1
    W=who;
    out_var(W{:});
end

%% %%%%%%%%%%%%%%%%%%% figure 8 %%%%%%%%%%%%%%%%%%%%%

f8=my_fig(8);
copyobj(get(f4.s1,'children'),f8.s1);
set(f8.s1,'xlim',[0 0.1],'ylim',[-0.08 0]);
center_axes(f8.s1);

%% %%%%%%%%%%%%%%%%%%% figure 9 %%%%%%%%%%%%%%%%%%%%%

f9.f=figure(9); clf(f9.f);
f9.s1=copyobj(findall(f6.f,'type','axes'),f9.f);
f9.s1=findall(f9.f,'type','axes');
delete(findall(f9.f,'type','surface'));
scatter3(f9.s1,IA(m3,7),IA(m3,6),m3_m,'filled','cdata',m3_m,'sizedata',5);


disp(['fini: ',datestr(clock)]);


function out_var(varargin)
% This function output the function variable space to the base workspace
for dum=1:numel(varargin)
    assignin('base',varargin{dum},evalin('caller',varargin{dum}));
end


function I=clean(I,ii)
% This function removes unwanted pixels in the image by designating nans.
% 'ii' represents the indices corresponding to pxs that will be replaced
% with nans
I1=I(:,:,1);%red
I2=I(:,:,2);%green
I3=I(:,:,3);%blue

I1(ii)=nan;
I2(ii)=nan;
I3(ii)=nan;
I=cat(3,I1,I2,I3);


function y=mypoly2(x,p)
y=nan(size(x));    
y(x<p(2))=0;
y(x>=p(2))=p(1).*(x(x>=p(2))-p(2)).^2;

