function parse_sp_crack02(filename,frame_n,varargin)
% Author: Joshua Yeh
% Date created: 2018-04-18
% Updated on: 2018-05-15
% 
%% DESCRIPTION FOR UNLOADING ANALYSIS
% This script accepts a video frame and analyzes the pixels in the frame.
% The pixels are categorized as loading, unloading, background, or other.
% The pixels are plotted on a TCC vs GCC chromatic plot.
% 
%% INPUT VARIABLES
% filename: the name of the file that will be used to perform the analysis
% 
% frame_n: the frame in which the analysis will be performed on
% 
% varargin: <'field_name'>,<'value'>
    % 'skip': A toggle flag where a value of 1 signifies to the function to
    % skip manual import of the frame from the video file. A value of 0
    % signifies the opposite.
    %
    % 'ref_frame': the reference frame representing the beginning of the
    % experiment (default is 1)
    %
    % 'ROI': the region of interest in the main frame
    % 
    % 'ref_ROI': the region go interest in the reference frame
    %
    % 'white_ROI': the region of interest representing the white area in
    % the main frame
    %
    % 'white_ref_ROI': the region of interest representing the white area
    % in the reference frame
    %
    % 'background_thresh': the background threshold value that
    % distinguishes the background from the foreground (for masking
    % purposes)
    %
    % 'dark_thresh': the dark thereshold value that excludes the dark
    % regions in the sample (for masking purposes)
    %
    % 'stress_calc': toggle (1 or 0, 1 is default) for performing stress
    % association of identified loading and unloading responses
    %
    % 'export_var': export variables from the function workspace to the
    % base workspace
    %
    % 'TCC_thresh': percentile (0 to 1) threshold distinguishing signal
    % from noise
    %
    % 'GCC_thresh': percentil (0 to 1) threshold distinguishing sigal from
    % noise
    %
    % 'import_mask': this is a flag (0 (default) or 1) to manually import
    % an image mask and use it in addition to the masking algorithm written
    % in this function
    %
    % 'import_mask_file': this is the filename (be sure to include the
    % path) of the image mask the user can import and use
    %
    %
    % 'noise_filter': this is a flag (0 (default) or 1) for applying a
    % noise filter after the picel parsing process
    %
    %
    % 'video_mode': this is a flag (0 (default) or 1) for running the
    % algorithm in video mode
    %
    % 'res': the step size for determining the x and y edges for the
    % hexagonal histogram binning and plotting for the chromatic density
    % plots
    % 
    % 'max': the maximum counts at saturation for the honeycomb 2d
    % histogram
    % 
    % 'px2mm': pixel to mm resolution
    %
    % 'MCvps': Merocyanine to peak stress calibration curve
    %
    % 'MC_calc': toggle for calculating MC map (0 is default)
    %
    % 'ED_calc': toggle for calculating the energy density (0 is default)
    %
    % 'ED_unload': mat file containing data that converts unloading stress
    % values to energy density values
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
params.addParameter('white_ROI',[600 700 350 450],...
    @(x) isnumeric(x)&numel(x)==4);
params.addParameter('white_ref_ROI',[800 1000 500 600],...
    @(x) isnumeric(x)&numel(x)==4);
params.addParameter('background_thresh',565,@(x) isnumeric(x));
params.addParameter('dark_thresh',340,@(x) isnumeric(x));
params.addParameter('stress_calc',1,@(x) islogical(x)|x==1|x==0);
params.addParameter('export_var',1,@(x) islogical(x)|x==1|x==0);
params.addParameter('TCC_thresh',0.95,@(x) isnumeric(x)&x>0&x<1);
params.addParameter('GCC_thresh',0.05,@(x) isnumeric(x)&x>0&x<1);
params.addParameter('import_mask',0,@(x) isnumeric(x));
params.addParameter('import_mask_file',' ',@(x) ischar(x));
params.addParameter('noise_filter',0,@(x) isnumeric(x));
params.addParameter('video_mode',0,@(x) isnumeric(x));
params.addParameter('initialize','normal',...
    @(x) ischar(x)&(strcmp(x,'only')|strcmp(x,'normal')));
params.addParameter('res',2e-3,@(x) isnumeric(x));
params.addParameter('max',10,@(x) isnumeric(x)&x>0);
params.addParameter('px2mm',38.94,@(x) isnumeric(x));
params.addParameter('MCvps','MCvps_TN135_180430.mat');
params.addParameter('MC_calc',0,@(x) isnumeric(x));
params.addParameter('ED_calc',0,@(x) isnumeric(x));
params.addParameter('ED_unload','');
params.parse(varargin{:});

%Extract out values from parsed input
skip=params.Results.skip;
ref_frame=params.Results.ref_frame;
ROI=params.Results.ROI;
ref_ROI=params.Results.ref_ROI;
white_ROI=params.Results.white_ROI;
white_ref_ROI=params.Results.white_ref_ROI;
background_thresh=params.Results.background_thresh;
dark_thresh=params.Results.dark_thresh;
stress_calc=params.Results.stress_calc;
export_var=params.Results.export_var;
import_mask=params.Results.import_mask;
noise_filter=params.Results.noise_filter;
video_mode=params.Results.video_mode;
init=params.Results.initialize;
res=params.Results.res;
max1=params.Results.max;
px2mm=params.Results.px2mm;
MCvps=params.Results.MCvps;
MC_calc=params.Results.MC_calc;
ED_calc=params.Results.ED_calc;
ED_unload=params.Results.ED_unload;

if video_mode==1
    visible='off';
else
    visible='on';
end

%% LOAD VIDEO AND FITS
disp(['fcn performing on frame: ',num2str(frame_n)]);
[~,name,~]=fileparts(filename);
if skip==0
    frame=extract_frames(filename,[ref_frame frame_n]);
    save(['mat_files/UL_',name,'_',num2str(frame_n),'.mat'],'frame');
else
    try
        load(['mat_files/UL_',name,'_',num2str(frame_n),'.mat']);
    catch
        disp('mat file was not found, attempting to extract from file');
        frame=extract_frames(filename,[ref_frame frame_n]);
        save(['mat_files/UL_',name,'_',num2str(frame_n),'.mat'],'frame');
    end
end


% Load the loading fit curves for TN135

%fits are based on data from: F:\Dropbox\ESPCI-Paris\My
%Research\data\Yinjun\Cyclic loading\20180430\TN135_replot.m
load('TN135_loading_fits.mat','TCCvGCC_load','temp_TCC_load_avg',...
    'temp_TCC_load_std','load_fit2_TN135');

disp('Frame loaded!');

if strcmp(init,'only')
    disp(['Fini: ',datestr(clock)]);
    return
end

%% PERFORM COLOR CORRECTIONS

% Crop our ROI
subimage=frame(2).CData(ROI(1):ROI(2),ROI(3):ROI(4),:);
ref=frame(1).CData(ref_ROI(1):ref_ROI(2),ref_ROI(3):ref_ROI(4),:);

% Define white refference area
white=frame(2).CData(white_ROI(1):white_ROI(2),white_ROI(3):white_ROI(4),:);
white_ref=frame(1).CData(white_ref_ROI(1):white_ref_ROI(2),...
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
ii=I1>background_thresh;
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
template(~isnan(frame_corr(:,:,1)))=1;%set elements that are #s as 1s
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
RCC=R_ratio-nanmean(R_ratio0(:));
GCC=G_ratio-nanmean(G_ratio0(:));
BCC=B_ratio-nanmean(B_ratio0(:));
TCC=sqrt(RCC.^2+GCC.^2+BCC.^2);

% reference frame
RCC0=R_ratio0-nanmean(R_ratio0(:));
GCC0=G_ratio0-nanmean(G_ratio0(:));
BCC0=B_ratio0-nanmean(B_ratio0(:));
TCC0=sqrt(RCC0.^2+GCC0.^2+BCC0.^2);


%% CATEGORIZE PIXELS AS UNLOADING, LOADING, OTHER, OR BACKGROUND

% import the 0-iso stress line from color-stress map for TN135
try delete(fcsm.f); catch; end
fcsm.f=openfig('CSM_TN135.fig','invisible');
fcsm.xdata=get(findall(fcsm.f,'userdata','zero_stress'),'xdata');
fcsm.ydata=get(findall(fcsm.f,'userdata','zero_stress'),'ydata');
try delete(fcsm.f); catch; end

% Create indexed array of del_G_ratio and TCC
disp('Parsing color coordinates');
TCC_n=(TCC(:)-temp_TCC_load_avg)./temp_TCC_load_std;%normalized TCC
load_col=polyval(TCCvGCC_load,TCC_n);
IA=[(1:numel(TCC(:)))',GCC(:),TCC(:),load_col,...
    zeros(numel(TCC(:)),1),...%flag col
    spline(fcsm.ydata,fcsm.xdata,TCC(:)),...% GCC of 0 iso stress line
    spline(fcsm.xdata,fcsm.ydata,GCC(:)),...% TCC of 0 iso stress line
    ones(numel(TCC(:)),1),ones(numel(TCC(:)),1)];

% Determine the closest distance between the current pixel chromatic
% coordinate and the chromatic loading curve.
[~,d]=knnsearch(IA(:,[4 3]),IA(:,2:3));

% Get threshold values
thresh_GCC=quantile(GCC0(:),params.Results.GCC_thresh);
thresh_TCC=quantile(TCC0(:),params.Results.TCC_thresh);

for dum=1:numel(TCC)
    [r,c]=ind2sub(size(TCC),dum);
    IA(dum,8:9)=[r,c];
    
    if isnan(IA(dum,2))
        continue
    end
        
    % This corresponds to the unloading regions.
    if d(dum)>abs(thresh_GCC) && IA(dum,2)<IA(dum,4) &&...
            IA(dum,2)>IA(dum,6) && IA(dum,3)>IA(dum,7)...
            &&IA(dum,3)>thresh_TCC
        IA(dum,5)=1;
        
    % This corresponds to the background noise.
    elseif IA(dum,2)>thresh_GCC&&IA(dum,3)<thresh_TCC
        IA(dum,5)=2;
    
    % This corresponds to the loading regions
    elseif d(dum)<abs(thresh_GCC) &&...
            IA(dum,3)>thresh_TCC
        IA(dum,5)=3;
        
    end
    if isnan(IA(dum,2))%don't flag nan values
        IA(dum,5)=nan;
    end
    
    if mod(dum,1e4)==0
        disp([num2str(dum),' of ',num2str(numel(TCC)),' parsed']);
    end
end
disp([num2str(dum),' of ',num2str(numel(TCC)),' parsed']);

% Extract the indices associated with each category
m3=find(IA(:,5)==3);%indices associated with loading color change
m2=find(IA(:,5)==2);%indices associated with sample background noise
m1=find(IA(:,5)==1);%indices associated with unloading color change
m0=find(IA(:,5)==0);%indices not associated with unloading/loading color change
nan1=find(isnan(IA(:,5)));%indices associated with background

% NOISE FILTERING
if noise_filter==1%if flag is turned on, apply a noise filter to the parsed pxs
    % Apply noise filter on a smaller local scale (unloading). This removes
    % isolated pixels.
    template=zeros(size(TCC));
    template(m1)=1;%set elements corresponding to unloading to 1
    local_sum=conv2(template,ones(3),'same');%compute local sum
    ii=local_sum<=2&local_sum>0;
    IA(ii,5)=0;
    % Refresh loading points
    m1=IA(:,5)==1;%indices associated with unloading color change

    % Apply noise filter on a smaller local scale (loading) This removes
    % isolated pixels.
    template=zeros(size(TCC));
    template(m3)=1;%set elements corresponding to unloading to 1
    local_sum=conv2(template,ones(3),'same');%compute local sum
    ii=local_sum<=2&local_sum>0;
    IA(ii,5)=0;
    % Refresh loading points
    m3=IA(:,5)==3;%indices associated with loading color change

    % Apply noise filter on a larger local scale
    template=zeros(size(TCC));
    template(m1)=1;%set elements corresponding to unloading to 1
    template(m3)=1;%set elements corresponding to unloading to 1
    local_sum=conv2(template,ones(8),'same');%compute local sum
    ii=local_sum<=5&local_sum>0;
    IA(ii,5)=0;
    % Refresh loading and unloading points
    m1=find(IA(:,5)==1);%indices associated with unloading color change
    m3=find(IA(:,5)==3);%indices associated with loading color change
end

%% ASSOCIATE UNLOADING AND LOADING POINTS WITH STRESS VALUES

if MC_calc==1
    MCvps0=load(MCvps);%load the MC vs stress calibration data
end

if stress_calc==1
    load('TN135_color_stress_map.mat','isrx','isry','stress_q2_TN135','urx',...
    'ury');
    isrx2=isrx(:);%reform to column array to speed up process
    isry2=isry(:);%reform to column array
    urx2=urx(:);% reform to col array
    ury2=ury(:);% reform to col array
    
    % Create loading curve array (stress vs. TCC)
    ss2=linspace(0,10,1e3)';%stress values associated with loading curve
    TTC2=load_fit2_TN135(ss2);%TCC values associated with loading curve
    
    disp('Associating unloading points to stresses');
    m1_s=nan(numel(m1),1);%variable containing stress associated with unloading
    m1_ps=m1_s;%variable containing the peak stress associated with unloading
    m1_mc=m1_s;%variable containing the [MC] associated with the peak loading
    m1_ed=m1_s;%variable containing the strain energy density of unload. px.
    
    % process unloading points
    for dum=1:numel(m1)
        del=abs(IA(m1(dum),2)-isrx2).^2+abs(IA(m1(dum),3)-isry2).^2;
        [~,c]=ind2sub(size(isrx),find(del==min(del),1));
        m1_s(dum)=stress_q2_TN135(c);%stress in MPa
        
        % Associate chromatic points with peak stress value (MPa) by
        % tracing the unloading curve back to the loading curve.
        del=abs(IA(m1(dum),2)-urx2).^2+abs(IA(m1(dum),3)-ury2).^2;
        [~,c]=ind2sub(size(urx),find(del==min(del),1));
        cx=urx(end,c);%GCC unloading values @ peak stress
        cy=ury(end,c);%TCC unloading values @ peak stress
        del=abs(cy-TTC2);
        ii=find(del==min(del),1);
        m1_ps(dum)=ss2(ii);%peak stress associated with unloading (MPa)
        
        % Convert the peak stress value to [MC] (%)
        if MC_calc==1
            m1_mc(dum)=mypoly2(m1_ps(dum),MCvps0.MCvps);%MC conc (%)
        end
        
        % Convert stress values to energy density values
        if ED_calc==1
            ED_data=load(ED_unload);%load the mat needed to perform calc.
            psq=m1_ps(dum);%peak stress associated with unloading (MPa)
            
            % det. unloading fit parameters for current peak stress value
            par_q1=nan(numel(fieldnames(ED_data.UPMF)),1);
            for dum2=1:numel(par_q1)
                par_q1(dum2)=ED_data.UPMF.(['order',num2str(dum2-1)])(psq);
            end
            
            
            % build interpolated unloading stress vs. \lambda curve
            unload_LAM=linspace(1,ED_data.load_fit2(psq),1000)';% lambda
            
            % create normalized lambda for fit based on normalized xdata
            unload_LAM_n=(unload_LAM-ED_data.f_Ux_mean(psq))./...
                ED_data.f_Ux_std(psq);
            %interp. unload. stress array
            unload_STRESS=polyval(par_q1,unload_LAM_n);
            
            
            % find the associated \lambda on the unloading curve based on
            % the unloading pixel's current stress
            del=abs(unload_STRESS-m1_s(dum));
            uu1=find(del==min(del));
            
            % Calculate strain energy density based on numerical integration
            m1_ed(dum)=trapz(unload_LAM(1:uu1),unload_STRESS(1:uu1));
            
        end

        if mod(dum,1e4)==0%update user on status
            disp([num2str(dum),' of ',num2str(numel(m1))]);
        end
    end
    disp('Finished associating unloading points to stress');

    disp('Associating loading points to stresses');
    m3_s=nan(numel(m3),1);%variable containing stress associated with loading
    m3_mc=m3_s;%variable containing MC concentration (calc from loading stress)
    m3_ed=m3_s;%variable containing the strain energy density from load. pxs.
    
    % If energy density will be calc., construct stress-lambda curves
    if ED_calc==1
        load_LAM=linspace(1,2.5,1e3);%lambda for loading curve
        load_STRESS=ED_data.load_fit(load_LAM);%loading stress (MPa)
    end
    
    % process loading points
    for dum=1:numel(m3)
        del=abs(IA(m3(dum),3)-TTC2);
        ii=find(del==min(del),1);
        m3_s(dum)=ss2(ii);%loading stress in MPa
        
        % Convert the peak stress value to [MC] (%)
        if MC_calc==1
            m3_mc(dum)=mypoly2(m3_s(dum),MCvps0.MCvps);%MC conc (%)
        end
        
        % Convert stress values to energy density values
        if ED_calc==1
            
            del=abs(load_STRESS-m3_s(dum));
            uu2=find(del==min(del));
            m3_ed(dum)=trapz(load_LAM(1:uu2),load_STRESS(1:uu2));%energy density
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

% Figure 1
f1=my_fig(1,{[1 2 1] [1 2 2]},'visible',visible);
axis([f1.s1 f1.s2],'image');
set(f1.f,'name','Raw extracted frame');

% Figure 2
f2=my_fig(2,{[3 2 1] [3 2 2] [3 2 3] [3 2 4] [3 2 5] [3 2 6]},...
    'fontsize',10,'visible',visible);
set(f2.f,'name','Corrected extracted frame');
title(f2.s1,'histogram of reference_{corr}','fontsize',12);
title(f2.s2,'histogram of frame_{corr}','fontsize',12);
title(f2.s3,'histogram of \DeltaG_{ratio} and \DeltaC_{tot} for ref',...
    'fontsize',12);
title(f2.s4,'histogram of \DeltaG_{ratio} and \DeltaC_{tot} for frame',...
    'fontsize',12);
title(f2.s5,'histogram of ref intensities');
title(f2.s6,'histogram of frame intensities');
set(f2.s3,'xlim',[-0.02 0.02]);
set(f2.s4,'xlim',[-0.02 0.02]);

% Figure 3
f3=my_fig(3,{[1 1 1]},'visible',visible);
f3.f.Position(4)=565;
f3.s2=axes;
f3.s3=axes;
set(f3.s2,'position',[0.65 0.7 0.1 0.3]);
set(f3.f,'name','Color mapping of reference sample area');
xylabels(f3.s1,'Green chromatic change','Total chromatic change');
set(f3.s3,'position',f3.s1.Position,'color','none','box','off',...
    'xcolor','none','ycolor','none','view',f1.s1.View);

% Figure 4
f4=my_fig(4,{[1 1 1]},'visible',visible);
f4.f.Alphamap=linspace(0,1,256);
set(f4.f,'name','Color mapping of activated sample area');
xylabels(f4.s1,'Green chromatic change','Total chromatic change');
set(f4.s1,'xlim',[-0.04 0.04],'ylim',[0 0.08]);

% Figure 5
f5=my_fig(5,{[1 1 1]},'marg_w',[0.01 0.01],'marg_h',[0.01 0.01],...
    'visible',visible);
f5.f.Alphamap=linspace(0,1,256);
f5.f.Position(3:4)=[667 645];

% Figure 6
f6=my_fig(6,{[1 2 1] [1 2 2]},'visible',visible);
f6.f.Alphamap=linspace(0,1,256);
f6.f.Position(2:4)=[377 959 463];
f6.s2.Position([2 4])=[0.4 0.55];
f6.s3=copyobj(f6.s2,f6.f);
uistack(f6.s3,'bottom');
xylabels(f6.s2,'Green chromatic change','Total chromatic change');

% Figure 7
if stress_calc==1
    f7=my_fig(7,{[1 2 1] [1 2 2]});
    f7.f.Position(3:4)=[959 660];
    f7.s2.Position([2 4])=[0.4 0.55];
    f7.s4=copyobj(f7.s1,f7.f);
    uistack(f7.s4,'top');
    xylabels(f7.s2,'Green chromatic change','Total chromatic change');
end

% Figure 8
f8=my_fig(8);
f8.f.Position(3:4)=[340 330];
axis(f8.s1,'image');

% Figure 9
f9=my_fig(9);
f9.f.Position(3:4)=[340 330];
f9.s2=copyobj(f9.s1,f9.f);

% Figure 10
f10=my_fig(10);

% Figure 11
f11=my_fig(11);
f11.f.Position(3:4)=[340 330];

% Figure 12
f12=my_fig(12);
f12.f.Position(3:4)=[340 330];

% Figure 13
f13=my_fig(13);
f13.f.Position(3:4)=[340 330];

% Figure 14
f14=my_fig(14);
f14.f.Position(3:4)=[340 330];

% import the iso-stress lines for TN135
try delete(fcsm.f); catch; end
fcsm.f=openfig('CSM_TN135.fig','invisible');
fcsm.s0=findall(fcsm.f,'type','axes');
clim_iso_stress=fcsm.s0.CLim;
delete(findall(fcsm.f,'color',ones(1,3).*0.7));
% Copy the plot to f3.f and f4.f figures
f3.csm=copyobj(findall(fcsm.s0(1),'type','line'),f3.s1);
f4.csm=copyobj(findall(fcsm.s0(1),'type','line'),f4.s1);
try delete(fcsm.f); catch; end

%% %%%%%%%%%%%%%%%%%%% figure 1 %%%%%%%%%%%%%%%%%%%%%
imagesc(f1.s1,uint8(ref_corr));
imagesc(f1.s2,uint8(frame_corr));

%% %%%%%%%%%%%%%%%%%%% figure 2 %%%%%%%%%%%%%%%%%%%%%

histogram(f2.s1,ref_corr(:,:,1),'facecolor','r','facealpha',0.5);
histogram(f2.s1,ref_corr(:,:,2),'facecolor','g','facealpha',0.5);
histogram(f2.s1,ref_corr(:,:,3),'facecolor','b','facealpha',0.5);
histogram(f2.s2,frame_corr2(:,:,1),'facecolor','r','facealpha',0.5);
histogram(f2.s2,frame_corr2(:,:,2),'facecolor','g','facealpha',0.5);
histogram(f2.s2,frame_corr2(:,:,3),'facecolor','b','facealpha',0.5);
histogram(f2.s3,TCC0(:),'facecolor','k','facealpha',0.5);
histogram(f2.s3,GCC0(:),'facecolor','g','facealpha',0.5);
histogram(f2.s4,TCC(:),'facecolor','k','facealpha',0.5);
histogram(f2.s4,GCC(:),'facecolor','g','facealpha',0.5);
histogram(f2.s5,ref_corr_sum(:));
histogram(f2.s6,frame_corr_sum(:));

%% %%%%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%
% Plot loading line
TCC_load=linspace(0,0.133,100);
TCC_load_n=(TCC_load-temp_TCC_load_avg)./temp_TCC_load_std;
plot3(f3.s1,polyval(TCCvGCC_load,TCC_load_n),TCC_load,...
    ones(1,numel(TCC_load)).*100,'k--');

% Create boundary of physical values on chromatic map
run G_RGB_mod_space;%plot boundary of color map
plot(f3.s1,G_space,RGB2,'-','color',ones(1,3).*0.7);
plot3(f3.s1,f3.s1.XLim,ones(1,2).*thresh_TCC,...
    ones(1,2).*100,'--','color',ones(1,3).*0.7);
plot3(f3.s1,polyval(TCCvGCC_load,TCC_load_n)+thresh_GCC,TCC_load,...
    ones(1,numel(TCC_load)).*100,'--','color',ones(1,3).*0.7);
plot3(f3.s1,polyval(TCCvGCC_load,TCC_load_n)-thresh_GCC,TCC_load,...
    ones(1,numel(TCC_load)).*100,'--','color',ones(1,3).*0.7);
f3.c1=colorbar(f3.s1,'location','eastoutside');
f3.c1.YLabel.String='Stress (MPa)';
colormap(f3.s1,'summer');
set(f3.s1,'clim',clim_iso_stress,'xlim',[-0.1 0.05],'ylim',[0 0.15]);
center_axes(f3.s1);


%transform chromatic datapoints to 2d histogram array
% f3.surf1=my_hist2(f3.s1,GCC0(:),TCC0(:));
f3.surf1=my_honeycomb(f3.s1,GCC0(:),TCC0(:),'max',100,'res',res,...
    'dn','background','color','k');

view(f3.s3,[0 90]);
xylabels(f3.s3,' ',' ');
set(f3.s3,'xlim',f3.s1.XLim,'ylim',f3.s1.YLim,'xtick',[],'ytick',[],...
    'clim',[0 100],'view',f3.s1.View,'position',f3.s1.Position,...
    'fontname',f3.s1.FontName,'fontsize',f3.s1.FontSize);
f3.c3=alphaColorbar(f3.s3,f3.surf1,'clab','Counts');
f3.s1.Position=f3.s3.Position;
linkaxes([f3.s1 f3.s3],'xy');


% image frame of reference area
imagesc(f3.s2,uint8(ref_corr));
axis(f3.s2,'image');
set(f3.s2,'box','off','xcolor','none','ycolor','none','xtick',[],...
    'ytick',[]);



%% %%%%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%
plot(f4.s1,polyval(TCCvGCC_load,TCC_load_n),TCC_load,'k--');

% Create boundary of physical values on chromatic map
plot(f4.s1,G_space,RGB2,'-','color',ones(1,3).*0.7);
copyobj(findall(f3.s1,'type','line','linestyle','--'),f4.s1);
set(f4.s1,'clim',clim_iso_stress,'ylim',[0 0.15],'xlim',[-0.1 0.05]);
f4.c1=colorbar(f4.s1,'location','eastoutside');
f4.c1.YLabel.String='Stress (MPa)';
colormap(f4.s1,'summer');
center_axes(f4.s1);

% Transform chromatic datapoints to 2d histogram array (m0)
f4.N_m0=my_honeycomb(f4.s1,GCC(m0),TCC(m0),'max',max1,'res',res,...
    'dn','unloading','color','k');
f4.N_m0.ZData=ones(size(f4.N_m0.XData)).*100;

% Transform chromatic datapoints to 2d histogram array (m1)
f4.N_m1=my_honeycomb(f4.s1,GCC(m1),TCC(m1),'max',max1,'res',res,...
    'dn','unloading','color','m');
f4.N_m1.ZData=ones(size(f4.N_m1.XData)).*100;

% Transform chromatic datapoints to 2d histogram array (m3)
f4.N_m3=my_honeycomb(f4.s1,GCC(m3),TCC(m3),'max',max1*1e1,'res',res,...
    'dn','loading','color','b');
f4.N_m3.ZData=ones(size(f4.N_m3.XData)).*100;


%% %%%%%%%%%%%%%%%%%%% figure 5 %%%%%%%%%%%%%%%%%%%%%
imagesc(f5.s1,uint8(rgb_correction(subimage,white,'simple',255)));
plot(f5.s1,IA(m3,9),IA(m3,8),'b.');
plot(f5.s1,IA(m1,9),IA(m1,8),'r.','tag','map');
patch(f5.s1,[0 78.94 78.94 0 0],[5 5 25 25 5],'w','edgecolor','none',...
    'facealpha',0.8);
plot(f5.s1,[0 48.94]+15,[10 10],'k-','linewidth',4);
text(f5.s1,32.25+10,10,'1 mm','horizontalalignment','center',...
    'verticalalignment','bottom');
axis(f5.s1,'image');
set(f5.s1,'box','off','xtick',[],'ytick',[],'xcolor','none',...
    'ycolor','none');


%% %%%%%%%%%%%%%%%%%%% figure 6 %%%%%%%%%%%%%%%%%%%%%

copyobj(get(f5.s1,'children'),f6.s1);% copy image frame
axis(f6.s1,'image');
f6.s1.Position(1)=0.07;
set(findall(f6.s1,'type','patch'),'ydata',[5 5 35 35 5]);
set(findall(f6.s1,'type','text'),'position',[40 10 0]);
set(f6.s1,'xtick',[],'ytick',[],'xcolor','none','ycolor','none')

copyobj(get(f4.s1,'children'),f6.s2);
set(findall(f6.s2,'type','line','linestyle','-'),'parent',f6.s3);
set(f6.s2,'xlim',[-0.04 0.04],'ylim',[0 0.08]);
linkaxes([f6.s2 f6.s3],'xy');

% colorbar
f6.c3=colorbar(f6.s3,'horizontal');
colormap(f6.s3,'summer');
set(f6.c3,'position',[0.6048 0.175 0.3452 0.0461]);
f6.c3.YLabel.String='Stress (MPa)';
set(f6.s3,'clim',clim_iso_stress,'position',f6.s2.Position,...
    'xlim',[-0.06 0.02]);

%% %%%%%%%%%%%%%%%%%%% figure 7 %%%%%%%%%%%%%%%%%%%%%

if stress_calc==1
    copyobj(findall(f6.s1,'type','image'),f7.s1);% copy image frame
    axis(f7.s1,'image');
    f7.s1.Position(1)=0.07;
    set(f7.s1,'xtick',[],'ytick',[],'xcolor','none','ycolor','none')
    
    % Plot chromatic coordinates on the image frame
    scatter3(f7.s4,IA(m1,9),IA(m1,8),m1_s,'filled','cdata',m1_s,...
        'sizedata',5,'markeredgecolor','none');
    colormap(f7.s4,spring);
    scatter3(f7.s1,IA(m3,9),IA(m3,8),m3_s,'filled','cdata',m3_s,...
        'sizedata',5,'markeredgecolor','none');
    colormap(f7.s1,flipud(winter));
    axis(f7.s4,'image');
    
    % find the stress associated at the TCC threshold
    set(f7.s1,'clim',[1 5]);
    set(f7.s4,'clim',[0 1.5]);
    
    % loading colorbar
    f7.c1=colorbar(f7.s1,'horizontal');
    f7.s1.Position([2 4])=[0.4 0.55];
    f7.c1.Position(1:4)=...
        [f7.s1.Position(1) 0.2455,...
        f7.s1.Position(3) 0.05];
    f7.c1.YLabel.String='virgin loaded';
    f7.c1.YLabel.Position(2)=0.9995;
    set(f7.c1,'box','off','tickdir','both','limits',[0 6]);
    
    % unloading colorbar
    f7.c4=colorbar(f7.s4,'horizontal');
    f7.c4.Position(1:4)=[f7.s1.Position(1) 0.1227,...
        f7.s1.Position(3) 0.05];
    set(f7.c4,'yaxislocation','top','yticklabel',[],'tickdir','both',...
        'box','off','limits',[0 6]);
    f7.c4.YLabel.String='unloaded';
    f7.c4.YLabel.Position(2)=0.0220;
    axes(f7.s4);
    f7.ct=annotation('textbox');
    f7.ct.Position(1:2)=[f7.s1.Position(1)+...
        (f7.s1.Position(3)-f7.ct.Position(3))/2 0.275];
    set(f7.ct,'horizontalalignment','center','string','Stress (MPa)',...
        'edgecolor','none','backgroundcolor','none',...
        'fontname','Microsoft YaHei Light');
    set(f7.s4,'position',f7.s1.Position,'xlim',f7.s1.XLim,'ylim',f7.s1.YLim,...
        'box','off','xcolor','none','ycolor','none','xtick',[],'ytick',[]);
    linkaxes([f7.s1 f7.s4],'xy');
    
    % chromatic-stress mapping
    f7.s3=copyobj(f6.s3,f7.f);
    copyobj(get(f4.s1,'children'),f7.s2);
    delete(findall(f7.s2,'type','line','linestyle','-'));
    set(f7.s2,'xlim',[-0.04 0.04],'ylim',[0 0.08]);
    uistack(f7.s3,'bottom');
    linkaxes([f7.s2 f7.s3],'xy');
    f7.c3=colorbar(f7.s3);
    colormap(f7.s3,'summer');
    f7.c3.YLabel.String='Stress (MPa)';
    set(f7.s3,'clim',clim_iso_stress,'position',f7.s2.Position,...
        'xlim',[-0.06 0.02]);
    f7.s3.Position(1)=0.55;
    f7.s3.Position(2)=0.45;
    f7.s3.Position(3)=0.3;
    f7.s3.Position(4)=0.5;
    f7.s2.Position=f7.s3.Position;
      
    % create custom colorbar axes for 2d histograms
    f7.s7=alphaColorbar_stacked(f7.s3,max1);
    xpos=10^(log10(f7.s7.XLim(2))/2);
    f7.t1=text(f7.s7,xpos,0.6,'other objects');
    set(f7.t1,'fontname',f7.s7.FontName,'horizontalalignment','center',...
        'color','w');
    f7.t2=text(f7.s7,xpos,1.6,'unloading');
    set(f7.t2,'fontname',f7.s7.FontName,'horizontalalignment','center');
    f7.t3=text(f7.s7,xpos,2.6,'loading');
    set(f7.t3,'fontname',f7.s7.FontName,'horizontalalignment','center');
    
%     linkaxes([f7.s3 f7.s5 f7.s6],'xy');
end

%% %%%%%%%%%%%%%%%%%%% figure 8 %%%%%%%%%%%%%%%%%%%%%
copyobj(findall(f6.s1,'type','image'),f8.s1);% copy image frame
set(f8.s1,'xtick',[],'ytick',[]);
text(f8.s1,7,20,'1 mm','fontname',f8.s1.FontName);
plot(f8.s1,[7 7+px2mm],[10 10],'k-','linewidth',2);
center_axes(f8.s1);
set(f8.s1,'layer','top');

%% %%%%%%%%%%%%%%%%%%% figure 9 %%%%%%%%%%%%%%%%%%%%%
% apparant, correlating uniaxial stress
copyobj(findall(f6.s1,'type','image'),f9.s1);% copy image frame
if stress_calc==1
    copyobj(findall(f7.s1,'type','scatter'),f9.s1);
    copyobj(findall(f7.s4,'type','scatter'),f9.s2);
    colormap(f9.s1,flipud(winter));
    colormap(f9.s2,spring);
end
axis(f9.s1,'image');
axis(f9.s2,'image');
set(f9.s1,'xtick',[],'ytick',[],'clim',f7.s1.CLim);
text(f9.s1,7,20,'1 mm','fontname',f8.s1.FontName);
plot(f9.s1,[7 7+px2mm],[10 10],'k-','linewidth',2);
center_axes(f9.s1);
set(f9.s1,'layer','top');
set(f9.s2,'xtick',[],'ytick',[],'clim',f7.s4.CLim,'xlim',f9.s1.XLim,...
    'ylim',f9.s1.YLim);
center_axes(f9.s2);
drawnow;
f9.s2.Position=f9.s1.Position;
f9.f.Name='Correlated uniaxial stress';

%% %%%%%%%%%%%%%%%%%%% figure 10 %%%%%%%%%%%%%%%%%%%%
plot(f10.s1,TTC2,ss2,'k-');
xylabels(f10.s1,'Total chromatic change','Stress (MPa)');
center_axes(f10.s1);

%% %%%%%%%%%%%%%%%%%%% figure 11 %%%%%%%%%%%%%%%%%%%%
% [MC] mapping
if MC_calc==1
    copyobj(findall(f6.s1,'type','image'),f11.s1);% copy image frame
    scatter3(f11.s1,IA(m1,9),IA(m1,8),m1_mc,'filled','cdata',m1_mc,...
        'sizedata',5,'markeredgecolor','none');
    scatter3(f11.s1,IA(m3,9),IA(m3,8),m3_mc,'filled','cdata',m3_mc,...
        'sizedata',5,'markeredgecolor','none');
    axis(f11.s1,'image');
    set(f11.s1,'xtick',[],'ytick',[],'layer','top','clim',[0 5]);
    center_axes(f11.s1);
end
f11.f.Name='[MC] (%)';

%% %%%%%%%%%%%%%%%%%%% figure 12 %%%%%%%%%%%%%%%%%%%%
if stress_calc==1
    % Amount of stress relaxed for unloading pixels
    copyobj(findall(f6.s1,'type','image'),f12.s1);% copy image frame
    scatter3(f12.s1,IA(m1,9),IA(m1,8),m1_ps-m1_s,'filled','cdata',m1_ps-m1_s,...
        'sizedata',5,'markeredgecolor','none');
    axis(f12.s1,'image');
    set(f12.s1,'xtick',[],'ytick',[],'layer','top','clim',[2 5]);
    f12.f.Name='Amount stress was relaxed (unloading pxs)';
    center_axes(f12.s1);
end

%% %%%%%%%%%%%%%%%%%%% figure 13 %%%%%%%%%%%%%%%%%%%%
if stress_calc==1
    % Peak loading stress for unloading pixels (loading history)
    copyobj(findall(f6.s1,'type','image'),f13.s1);% copy image frame
    scatter3(f13.s1,IA(m1,9),IA(m1,8),m1_ps,'filled','cdata',m1_ps,...
        'sizedata',5,'markeredgecolor','none');
    axis(f13.s1,'image');
    set(f13.s1,'xtick',[],'ytick',[],'layer','top','clim',[4 7]);
    f13.f.Name='Peak stress of unloading pixels';
    center_axes(f13.s1);
end

%% %%%%%%%%%%%%%%%%%%% figure 14 %%%%%%%%%%%%%%%%%%%%
% Energy density plots
if ED_calc==1&&stress_calc==1
    copyobj(findall(f6.s1,'type','image'),f14.s1);% copy image frame
    scatter3(f14.s1,IA(m1,9),IA(m1,8),m1_ed,'filled','cdata',m1_ed,...
        'sizedata',5,'markeredgecolor','none');
    scatter3(f14.s1,IA(m3,9),IA(m3,8),m3_ed,'filled','cdata',m3_ed,...
        'sizedata',5,'markeredgecolor','none');
    axis(f14.s1,'image');
    set(f14.s1,'xtick',[],'ytick',[],'layer','top','clim',[0 2]);
    f14.f.Name='Energy density mapping';
    center_axes(f14.s1);
end

%% Export function workspace to base workspace
if export_var==1
    W=who;
    out_var(W{:});
end

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
