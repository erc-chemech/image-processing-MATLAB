function f1=Argolight_IH(I1,test_I,varargin)
% Author: Joshua Yeh
% Date created: 2018-01-31
% 
%% DESCRIPTION
% IMPORTANT: THIS SCRIPT ONLY WORKS FOR IMAGES EXPORTED FROM THE NIKON
% AZ100 CONFOCAL MICROSCOPE.
% 
% This script performs an illumination inhomogeneity analysis on the
% ArgoLight field of rings (grid) pattern. This function only works on
% images exported by the Nikon AZ100 microscope.
% 
%% INPUT VARIABLES
% I1: structure variable that is exported using the import_tiff_stacks
% function
% 
% test_I: the image array in which the analysis will be based on
% 
% NAME PAIR ARGUMENTS: Argolight_gridfit(...'<fieldname>',<value>)
% 'T': threshold value for binarizing image (for identification of grid
% points) (default is T=100)
% 
%% OUTPUT VARIABLES
% 
%%

% Parse input variables
narginchk(1,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('T',100,@(x) isnumeric(x));
params.addParameter('span',0.05,@(x) isnumeric(x));
params.addParameter('pt_window',4,@(x) isnumeric(x));
params.addParameter('clim',[0 1],@(x) isnumeric(x)==1&numel(x)==2);
params.addParameter('dot_I','max',@(x) ischar(x)==1);
params.parse(varargin{:});

% Defines the subimage size around a grid pt. This helps reduce instances
% of multiple detected grid points associated with a single theoretical
% grid point placement.
pt_window=params.Results.pt_window;
clim=params.Results.clim;
dot_I=params.Results.dot_I;

% create anonymous fcns
% xy rotation matrix (counterclockwise) in degs
rotxy=@(x) [cosd(x) -sind(x);sind(x) cosd(x)];

% Binarize the average image
I1.bw=imbinarize(test_I,params.Results.T);
I1.bw=bwmorph(I1.bw,'fill');

% Extract the pixel resolution (um/px)
px_res=I1.info(1).UnknownTags(2).Value;
disp(['res: ',num2str(px_res),' ',956,'m/px']);

% Find location of grid points
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

% Construct grid assuming no rotation
spacing=50/px_res;%spacing between dots in pixels
R_row=centerA(1):spacing:I1.info(1).Width;
L_row=centerA(1):-spacing:1;
row=unique(sort([L_row R_row],'ascend'));
L_col=centerA(2):spacing:I1.info(1).Height;
U_col=centerA(2):-spacing:1;
col=unique(sort([U_col L_col],'ascend'));

AL_grid=nan(length(row)*length(col),2);
AL_grid2=nan(length(row)*length(col),2);%recenter grid at origin
count=1;
for dum1=1:length(row)
    for dum2=1:length(col)
        AL_grid(count,:)=[row(dum1) col(dum2)];
        AL_grid2(count,:)=[row(dum1)-centerA(1) col(dum2)-centerA(2)];
        count=count+1;
    end
end

% Find intensities corresponding to calibration dots
ii1=[];
for dum=1:size(centroids,1)
    c=[floor(centroids(dum,1)-pt_window) ceil(centroids(dum,1)+pt_window)];
    r=[floor(centroids(dum,2)-pt_window) ceil(centroids(dum,2)+pt_window)];
    if c(1)<1;   c(1)=1;     end
    if r(1)<1;   r(1)=1;     end
    if c(2)>size(test_I,2); c(2)=size(test_I,2);   end
    if r(2)>size(test_I,1); r(2)=size(test_I,1);   end
    subimage=test_I(r(1):r(2),c(1):c(2));
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

disp(['Constructed grid spacing: ',num2str(spacing)]);
disp(['Constructed ',num2str(size(AL_grid2,1)),' grid points']);
disp(['Detected ',num2str(size(ii1(:,1),1)+size(cross,1)),...
    ' grid points']);
disp(['Center grid pt is offset (rel. to theoretical image) center by: ',...
    newline,num2str(size(test_I)./2-centerA)]);

%% Fit a surface to the scatter datapoints

% Use this to perform loess surface fitting (local weighted quadratic
% fitting)
disp('Performing surface fit of intensity 3d scatter data');
ft = fittype( 'lowess' );
opts = fitoptions( 'Method', 'LowessFit',...
    'normalize','on','robust','bisquare','span',params.Results.span);

[sf, gof] = fit( [ii1(:,1), ii1(:,2)], ii1(:,3)./max(ii1(:,3)), ft, opts );
[inhomo_x,inhomo_y]=meshgrid(1:size(I1.tiff_stack,2),...
    1:size(I1.tiff_stack,1));

%% Plot data

% Create figure and format it
f1=my_fig(99,{[2 2 1] [2 2 2] [2 2 3] [2 2 4]},...
    'fontsize',8,'gap',[0.07 0.07],'marg_w',[0.1 0.01],'marg_h',[0.07 0.04]);
axis([f1.s1 f1.s2 f1.s4],'image');

set(f1.f,'units','normalized');
if f1.f.Position(2)>0.25
    f1.f.Position(2:4)=[0.1 0.5 0.7];
else
    f1.f.Position(3:4)=[0.5 0.7];
end


% subplot 1
surf(f1.s1,test_I);

% subplot 2
surf(f1.s2,test_I);
plot(f1.s2,ii1(:,1),ii1(:,2),'rx','displayname','grid pattern');
plot(f1.s2,centerA(1),centerA(2),'ro');
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
plot(sf);

%formatting
xylabels(f1.s1,'x','y');
xylabels(f1.s2,'x','y');
xylabels(f1.s3,'Intensity','Counts');
xylabels(f1.s4,'x','y');
colormap(f1.s1,'gray');
colormap(f1.s2,'gray');
set(f1.s4,'dataaspectratio',[1 1 0.0001],'clim',clim);
set([f1.s1 f1.s2 f1.s4],'ydir','reverse');
set(f1.s3,'yscale','log');
linkaxes([f1.s1 f1.s2 f1.s4],'xy');

% Output relevant data into f1
f1.sf=sf;%surface fit object
f1.center_offset=size(test_I)./2-centerA;% center offset displacement
f1.detected=size(ii1(:,1),1)+size(cross,1);% # grid pts detected
f1.constructed=size(AL_grid2,1);% # grid pts constructed
disp('Illumination inhomogeneity analysis completed');
