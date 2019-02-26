function f1=Argolight(I1,test_I,varargin)
% Author: Joshua Yeh
% Date created: 2018-02-02
% 
%% DESCRIPTION
% IMPORTANT: THIS SCRIPT ONLY WORKS FOR IMAGES EXPORTED FROM THE NIKON
% AZ100 CONFOCAL MICROSCOPE.
% 
% This script performs a illumination imhomogeneity and distortion analysis
% of a tif image exported by the Nikon AZ100 confocal microscope.
% 
%% INPUT VARIABLES
% I1: structure variable that is exported using the import_tiff_stacks
% function
% 
% test_I: the image array in which the analysis will be based on
% 
% 
% NAME PAIR ARGUMENTS: Argolight_gridfit(...'<fieldname>',<value>)
% 
% 'T': threshold value for binarizing image (for identification of grid
% points) (default is T=100)
% 
% 'span': span used for lowess or loess surface fitting procedure for
% inhomogeneous illumination
% 
% 'pt_window': the pixel distance from the centroid of a detected grid
% point, it is used to determine the appropriate itensity value assigned to
% the centroid of the detected grid point
%
%  'dot_I': the type of intensity assignment algorithm for grid point
%  detection, there are two types:
        % 'max' assigns the maximum intensity detected in the subimage to
        % the grid point
        % 'averaged' assigns the averaged intensities in the subimage with
        % intensity levels > 10 ADU (analogue to digital units)
% 
% 'clim': the color map span for the inhomo. illumination surface plot
% 
% 'fit_I': the type of fitting for inhomo. illum. surfac fitting
        % 'lowess' applies a quadratic loess fitting (default)
        % 'loess' applies a linear loess fitting
        % 'bihanrmonicinterp' applies a biharmonic surface interpolation
        % 'thinplateinterp' applies a thin plate interpolation
%         
% 'filter_outlier': option to remove outliers, as detected by isoutlier fcn
        % true (default)
        % false
% 'fig_n': figure number (default is 1) in which to plot the results in
%
%% OUTPUT VARIABLES
% f1: structure variables containing the foloowing fieldnames:
% 
    % f1.optimal_rot=optimal_rot;  fitted rotation of constructed grid
    % f1.AL_grid2=AL_grid2; coordinates of the constructed gridpoints
    % f1.centroids=centroids;   coordinates of the detected grid points
    % f1.cross=cross;   coordinates of the cross pattern gridpoint
    % f1.centerA=centerA;  coordinates of gridpoint closest to the image center
    % f1.sfcf=sfcf; linear fit for the conversion factor (photons to intensity)
    % f1.sfdx=sfdx; surface fit for distortion in the x direction
    % f1.sfdy=sfdy; surface fit for distortionin the y direction
    % f1.sfIH=sfIH; surface fit in illumination inhomogeneity
    % f1.I1=I1; image tiff structure variables that was used as the input
        % f1.I1.stats; information obtained from the region property fcn
        % applied on the binarized image of the test_I variable
    % f1.test_I=test_I; the image in which the analysis was performed on
    % f1.ii1; intensity for each corresponding detected gridpoint
    % f1.f; figure handles in which the results were plotted in
    % f1.s1; subplot 1 handle
    % f1.s2; subplot 2 handle
    % f1.s3; subplot 3 handle
    % f1.s4; subplot 4 handle
    % f1.s5; subplot 5 handle
    % f1.s6; subplot 6 handle
    % f1.s7; subplot 7 handle
    % f1.s8; subplot 8 handle
%     

%% Parse input variables

disp(' ');

narginchk(1,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('T',100,@(x) isnumeric(x));
params.addParameter('span',0.05,@(x) isnumeric(x));
params.addParameter('pt_window',4,@(x) isnumeric(x));
params.addParameter('clim',[0 1],@(x) isnumeric(x)==1&numel(x)==2);
params.addParameter('dot_I','averaged',@(x) ischar(x)==1);
params.addParameter('image_gain',2,@(x) isnumeric(x)&numel(x)==1);
params.addParameter('fit_I','lowess',@(x) ischar(x)==1);
params.addParameter('filter_outlier',true,@(X) islogical(x));
params.addParameter('fig_n',1,@(x) isnumeric(x));
params.parse(varargin{:});

% Extract out values from parsed input
pt_window=params.Results.pt_window;
clim=params.Results.clim;
dot_I=params.Results.dot_I;
image_gain=params.Results.image_gain;
fit_I=params.Results.fit_I;
filter_outlier=params.Results.filter_outlier;
fig_n=params.Results.fig_n;

% extract offset associated with image gain (defaut is 2)
% Dark offset values obtained from dark images. Indices corresponds to
% indices in the gain variable. These values were collected on 2018-02-01
offset=[6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 8 8 9 9 10];
C=offset(image_gain);

% Remove the offset from the image being used for analysis
test_I=test_I-C;

%% Create anonymous fcns
% xy rotation matrix (counterclockwise) in degs
rotxy=@(x) [cosd(x) -sind(x);sind(x) cosd(x)];

%% Binarize the average image
I1.bw=imbinarize(test_I,params.Results.T);
I1.bw=bwmorph(I1.bw,'fill');

%% Extract the pixel resolution (um/px)
px_res=I1.info(1).UnknownTags(2).Value;
disp(['res: ',num2str(px_res),' ',956,'m/px']);

%% Detect grid points
I1.stats=regionprops(I1.bw,'centroid','area','image','PixelList');
centroids0=reshape([I1.stats.Centroid],[2,size(I1.stats,1)])';

% Find the cross point (if any) by identifying outliers in the grid point
% area
TF=isoutlier([I1.stats(:).Area]);
iio=find(TF==1);
if ~isempty(iio)
    disp('identified cross-shaped point');
    out_area=[I1.stats(iio).Area];
    [~,iio2]=max(out_area);
    I1.stats2=I1.stats;
    I1.stats2(iio(iio2))=[];
    cross=[I1.stats(iio(iio2)).Centroid];
    centroids=reshape([I1.stats2.Centroid],[2,size(I1.stats2,1)])';
else
    disp('no cross-shaped point detected');
    cross=[nan nan];
    centroids=centroids0;
end

% Identify an approximate middle grid point
midx=size(test_I,2)/2;
midy=size(test_I,1)/2;
D=sqrt((centroids0(:,1)-midx).^2+(centroids0(:,2)-midy).^2);
[~,m3]=min(D);
centerA=centroids0(m3,:);

%% Construct grid assuming no rotation
spacing=50/px_res;%spacing between dots in pixels
R_row=centerA(1):spacing:I1.info(1).Width;
L_row=centerA(1):-spacing:1;
row=unique(sort([L_row R_row],'ascend'));
L_col=centerA(2):spacing:I1.info(1).Height;
U_col=centerA(2):-spacing:1;
col=unique(sort([U_col L_col],'ascend'));

AL_grid=nan(length(row)*length(col),2);
AL_grid2=nan(length(row)*length(col),2);
count=1;
for dum1=1:length(row)
    for dum2=1:length(col)
        %grid construction
        AL_grid(count,:)=[row(dum1) col(dum2)];
        
        %recenter constructed grid at origin
        AL_grid2(count,:)=[row(dum1)-centerA(1) col(dum2)-centerA(2)];
        count=count+1;
    end
end

%% Perform optimal rotation calculations
rots=linspace(-10,10,500);
e_rots=[];
for rot=rots
    % locate corresponding points
    pt_store=[];
    for dum=1:size(AL_grid,1)
        AL_grid_rot(dum,:)=(AL_grid(dum,:)-centerA)*rotxy(rot)+centerA;

        % find nearest point
        d=centroids-AL_grid_rot(dum,:);
        dtot=sqrt(d(:,1).^2+d(:,2).^2);
        [~,ii]=min(dtot);
        
        % determine the distance from the centerA reference point
        dc_total=norm(centroids(ii)-centerA);
        
        % store the calculated distances and weighted distances
        pt_store=[pt_store;ii dtot(ii) dtot(ii)*dc_total];
    end
    e_rots=[e_rots;nansum(pt_store(:,2)) rot];%store total displacement
end
[~,ii]=min(e_rots(:,1));%rotation w/ least total displacement is optimal fit
optimal_rot=rots(ii);
disp(['Optimal rotation (deg.): ',num2str(optimal_rot)]);
disp(['Constructed ',num2str(size(AL_grid2,1)),' grid points']);
disp(['Detected ',num2str(size(centroids(:,1),1)),...
    ' grid points']);
disp(['Center grid pt is offset (rel. to theoretical image) center by: ',...
    newline,num2str(size(test_I)./2-centerA)]);

% Rotate grid by the determined optimal rotation and recenter based on
% reference center location
for dum=1:size(AL_grid2,1)
    AL_grid2(dum,:)=AL_grid2(dum,:)*rotxy(optimal_rot)+centerA;
end

%% Distortion calculations
% Determine the x and y displacement of the reference rotated grid and the
% imaged grid pattern
DX=nan(size(AL_grid2,1),1);
DY=DX;
for dum=1:size(AL_grid2,1)
    dx=centroids(:,1)-AL_grid2(dum,1);%x-displacement
    dy=centroids(:,2)-AL_grid2(dum,2);%y-displacement
    dtot=sqrt(dx.^2+dy.^2);%total displacement
    
    % find the nearest point
    [~,ii]=min(dtot);
    
    if dtot(ii)<.7*spacing
        DX(dum)=dx(ii);
        DY(dum)=dy(ii);
    end
end

%% Find appropriate intensities corresponding to calibration dots
ii1=[];
for dum=1:size(centroids,1)
    
    % Define the bounds of the subimage
    c=[floor(centroids(dum,1)-pt_window) ceil(centroids(dum,1)+pt_window)];
    r=[floor(centroids(dum,2)-pt_window) ceil(centroids(dum,2)+pt_window)];
    if c(1)<1;   c(1)=1;     end
    if r(1)<1;   r(1)=1;     end
    if c(2)>size(test_I,2); c(2)=size(test_I,2);   end
    if r(2)>size(test_I,1); r(2)=size(test_I,1);   end
    
    subimage=test_I(r(1):r(2),c(1):c(2));%create the subimage
    [M,jj]=max(subimage(:));%find relative location of highest intensity
    [r1,c1]=ind2sub(size(subimage),jj);%find subscripts of high intensity
    
    kk0=find(subimage(:)>10);%find points associated with grid point
    averaged=mean(subimage(kk0));%find the mean intensity of the grid point
    
    % Determine absolute scaled magnitude distance from the image origin.
    % This is to provide a unique ID number associated with each grid
    % point.
    D=sqrt((c1+c(1))^2+(r1+r(1))^2);
    
    % Make sure the point is not associated with the cross-patterned point
    % by calculating the distance of point from the cross point
    if ~isnan(cross(1))
        DCA=norm([c1+c(1),r1+r(1)]-cross);
    else
        DCA=nan;
    end
    
    if ~isempty(ii1)
        kk=find(ii1(:,4)==D);
    else
        kk=[];
    end
    
    if isempty(kk)&&DCA>pt_window%only store point if it is unique
        if strcmp(dot_I,'max')
            ii1=[ii1;c1+c(1),r1+r(1),M,D];
        elseif strcmp(dot_I,'averaged')
            ii1=[ii1;centroids(dum,:),averaged,D];
        end
    end
end

if isempty(ii1)
    error('no points detected! analysis aborted');
end

%% Determine conversion factor (CF) between Intensity (ADU) and photons
disp('Performing photon conversion factor (CF) calculation');
n1=4;%super pixel size
n2=n1-1;
ydir=1:n1:size(test_I,1)-n1;% y-intervals defining super pixel
xdir=1:n1:size(test_I,2)-n1;% x-intervals defining super pixel

% Preallcoate large array
v_store=nan(numel(ydir)*numel(xdir),1);
m_store=v_store;

count=1;
flag=0;
for dum1=ydir
    for dum2=xdir
        %extract super pixel
        image_D=test_I(dum1:dum1+n2,dum2:dum2+n2);
        
        % Calculate variance of super pixel
        v_store(count)=var(image_D(:));
        m_store(count)=mean(image_D(:))-C;
        
        count=count+1;
    end
    ratio=dum1/size(test_I,1);
    if round(mod(ratio*100,25))==1&&round(mod(ratio*100,6))~=flag
        disp([num2str(ratio*100),'% completed']);
        flag=round(mod(ratio*100,6));
    end
end

disp('CF analysis complete');
%% Fit a surface to the scatter datapoints for inhomo. illum.

% Use this to perform loess surface fitting (local weighted quadratic
% fitting)
disp('Performing surface fit of intensity 3d scatter data');

switch fit_I
    case 'lowess'
        ft = fittype( 'lowess' );
        opts = fitoptions( 'Method', 'LowessFit',...
            'normalize','on','robust','bisquare','span',params.Results.span);
    case 'loess'
        ft = fittype( 'loess' );
        opts = fitoptions( 'Method', 'LowessFit',...
            'normalize','on','robust','bisquare','span',params.Results.span);
    case 'biharmonicinterp'
        ft = fittype( 'biharmonicinterp' );
        opts = fitoptions( 'Method', 'biharmonicinterp',...
            'normalize','on');
    case 'thinplateinterp'
        ft = fittype( 'biharmonicinterp' );
        opts = fitoptions( 'Method', 'biharmonicinterp',...
            'normalize','on');
end

% exclude outlier
if filter_outlier==true
    eo=isoutlier(ii1(:,3)./max(ii1(:,3)));
    eo=find(eo==0);
end

[sfIH, gof] = fit( [ii1(eo,1), ii1(eo,2)], ii1(eo,3)./max(ii1(eo,3)), ft, opts );
[inhomo_x,inhomo_y]=meshgrid(1:size(I1.tiff_stack,2),...
    1:size(I1.tiff_stack,1));

%% Prepare fitting procedure for distortion
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit','normalize','on','robust','bisquare');

% only use non-nan elements
ii=find(~isnan(DX));

% Fit a surface to x displacement field
disp('Performing x displacement surface fit');
[sfdx, gof] = fit( [AL_grid2(ii,1), AL_grid2(ii,2)], DX(ii), ft, opts );

% Fit a surface to y displacement field
disp('Performing y displacement surface fit');
[sfdy, gof] = fit( [AL_grid2(ii,1), AL_grid2(ii,2)], DY(ii), ft, opts );

%% Prepare fitting procedure for CF linear line
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares');
[sfcf, gof] = fit( m_store, v_store,ft, opts );
cf_coeff=coeffvalues(sfcf);

%% PLOT THE RESULTS
disp('Plotting the results...');

f1=my_fig(fig_n,...
    {[2 3 1] [2 3 2] [4 6 5:6] [2 3 4] [2 3 5] [4 6 17 23] [4 6 18 24],...
    [4 6 11:12]},...
    'fontsize',8,'gap',[0.12 0.1],'marg_w',[0.1 0.03],'marg_h',[0.09 0.04]);
axis([f1.s1 f1.s2 f1.s4 f1.s5 f1.s6 f1.s7],'image');

% subplot 1
imagesc(f1.s1,test_I);

% subplot 2
imagesc(f1.s2,test_I);
plot(f1.s2,ii1(:,1),ii1(:,2),'rx','displayname','grid pattern');
plot(f1.s2,centerA(1),centerA(2),'ro');
plot(f1.s2,AL_grid2(:,1),AL_grid2(:,2),'go');%show the rotated grid
plot(f1.s2,size(test_I,2)/2,size(test_I,1)/2,'m.','markersize',16);
if sum(isnan(cross))==0
    plot(f1.s2,cross(1),cross(2),'g+');
end

% subplot 3
histogram(f1.s3,test_I(:),logspace(0,3,1000));

% subplot 4
scatter3(f1.s4,ii1(:,1),ii1(:,2),ii1(:,3)./max(ii1(:,3)),...
    'cdata',ii1(:,3)./max(ii1(:,3)),...
    'markerfacecolor','flat','sizedata',10);
axes(f1.s4);
plot(sfIH);
sfIH_s=findall(f1.s4,'type','surface');
[~,sfIH_c3]=contour3(f1.s4,sfIH_s.CData,'showtext','on','linecolor','k');
set(sfIH_c3,'xdata',sfIH_s.XData,'ydata',sfIH_s.YData,'zdata',sfIH_s.ZData);

% subplot 5
imagesc(f1.s5,test_I);
quiver(f1.s5,AL_grid2(:,1),AL_grid2(:,2),DX,DY,1,'color','y',...
    'maxheadsize',0.5);

% subplot 6
scatter3(f1.s6,AL_grid2(:,1),AL_grid2(:,2),DX,'filled','cdata',DX);
axes(f1.s6);
plot(sfdx);

%subplot 7
scatter3(f1.s7,AL_grid2(:,1),AL_grid2(:,2),DY,'filled','cdata',DY);
axes(f1.s7);
plot(sfdy);

%subplot 8
plot(f1.s8,m_store,v_store,'kx');
axes(f1.s8);
plot(sfcf);
Th=text(f1.s8,0,0,['y=',num2str(cf_coeff(1)),'x+',num2str(cf_coeff(2))]);

% Axes formatting
set(Th,'units','normalized','position',[0.05 0.9 0],'fontsize',8)
set(f1.f,'name',I1.file);
set([f1.s1 f1.s2 f1.s4 f1.s5 f1.s6 f1.s7],'ydir','reverse');
f1.s8.XLim(1)=0;
f1.s8.YLim(1)=0;
title(f1.s1,'Original image');
title(f1.s2,'Detected and constructed FoR');
title(f1.s3,'Histogram of intensities');
title(f1.s4,'Illumination inhomogeneity');
title(f1.s5,'Distortion field');
title(f1.s6,'X-distortion');
title(f1.s7,'Y-distortion');
xylabels(f1.s1,'x','y');
xylabels(f1.s2,'x','y');
xylabels(f1.s3,'Intensity','Counts');
xylabels(f1.s4,'x','y');
xylabels(f1.s5,'x','y');
xylabels(f1.s6,'x','y');
xylabels(f1.s7,'x','y');
xylabels(f1.s8,'<I-offset>','var(I)');
colormap(f1.s1,'gray');
colormap(f1.s2,'gray');
colormap(f1.s5,'gray');
set(f1.s3,'yscale','log');
set(f1.s4,'dataaspectratio',[1 1 0.001]);
set(f1.s7,'dataaspectratio',[1 1 0.05]);
set(f1.s6,'dataaspectratio',[1 1 0.05]);
linkaxes([f1.s1 f1.s2],'xy');
drawnow;
disp('Plotting completed');

%% Store relevant information in output variables

f1.optimal_rot=optimal_rot;
f1.AL_grid2=AL_grid2;
f1.centroids=centroids;
f1.cross=cross;
f1.centerA=centerA;
f1.sfcf=sfcf;
f1.sfdx=sfdx;
f1.sfdy=sfdy;
f1.sfIH=sfIH;
f1.I1=I1;
f1.test_I=test_I;
f1.ii1=ii1;

disp('Argolight analysis completed');