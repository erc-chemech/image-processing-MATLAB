function f1=Argolight_gridfit(I1,test_I,varargin)
% Author: Joshua Yeh
% Date created: 2018/01/26
%% DESCRIPTION
% This script is only for tiff images taken with Nikon AZ100! Also, it is
% important that the calibration image is not severely rotated (no more
% than +/- 10 degrees in rotation). Also, for accurate results, be sure a
% grid point exist at the center of the image. Since distortion
% calculations are relative to the theoretical grid pattern (which is
% constructed from a center point extracted from the image), the
% theoretical grid pattern needs to be centered.
% 
% This script attempts to fit an Argolight field of rings pattern
% (Argorlight LM) to a calibration image This script also assumes that the
% tiff image that was imported was exported by the Nikon NIS Viewer. This
% is important in order to obtain accurate Metadata information. Grid
% points are spaced 50 microns apart.
% 
%% INPUT VARIABLES
% I1: structure variable that is exported using the import_tiff_stacks
% function
% 
% test: the image array in which the analysis will be based on
% 
% NAME PAIR ARGUMENTS: Argolight_gridfit(...'<fieldname>',<value>)
% 'T': threshold value for binarizing image (for identification of grid
% points) (default is T=100)
%     
%% OUTPUT VARIABLES
% f2: a structure variable containing fieldnames:
    % f2.s1: subplot 1 (upper left)
    % f2.s2: subplot 2 (upper left)
    % f2.s3: subplot 3 (upper right)
    % f2.s4: subplot 4 (upper right)
    % f2.s5: subplot 5 (lower left)
    % f2.s6: subplot 6 (lower right)
    % f2.I1: updated version of the I1 input argument
    % f2.optimal_rot: rotation applied to grid pattern
    % f2.DX: distortion of the image in the x direction
    % f2.DY: distortion of the image in the y direction
    % f2.AL_grid2: coordinates to the rotated grid pattern
    % f2.sfx: surface fit ofject for the distortion in the x direction
    % f2.sfy: surface fit object for the distortion in the y direction
%     
%%

% Parse input variables
narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('T',100,@(x) isnumeric(x));
params.parse(varargin{:});

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
centroids=reshape([I1.stats.Centroid],[2,size(I1.stats,1)])';

% Identify an approximate middle grid point
midx=size(test_I,2)/2;
midy=size(test_I,1)/2;
D=sqrt((centroids(:,1)-midx).^2+(centroids(:,2)-midy).^2);
[~,m3]=min(D);
centerA=centroids(m3,:);

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

% Perform optimal rotation calculations
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

% Determine the x and y displacement of the reference rotated grid and the
% imaged grid pattern
DX=nan(size(AL_grid2,1),1);
DY=DX;
for dum=1:size(AL_grid2,1)
    dx=centroids(:,1)-AL_grid2(dum,1);
    dy=centroids(:,2)-AL_grid2(dum,2);
    dtot=sqrt(dx.^2+dy.^2);%total displacement
    
    % find the nearest point
    [~,ii]=min(dtot);
    
    if dtot(ii)<.7*spacing
        DX(dum)=dx(ii);
        DY(dum)=dy(ii);
    end
end

% Prepare fitting procedure
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit','normalize','on','robust','bisquare');
% ft = fittype( 'poly55' );
% opts = fitoptions( 'Method', 'LinearLeastSquares','normalize','on',...
%     'robust','bisquare');

% only use non-nan elements
ii=find(~isnan(DX));

% Fit a surface to x displacement field
disp('Performing x displacement surface fit');
[sfx, gof] = fit( [AL_grid2(ii,1), AL_grid2(ii,2)], DX(ii), ft, opts );

% Fit a surface to y displacement field
disp('Performing y displacement surface fit');
[sfy, gof] = fit( [AL_grid2(ii,1), AL_grid2(ii,2)], DY(ii), ft, opts );


%% PLOT RESULTS
disp('Plotting results...');
f1=my_fig(99,{[2 2 1] [2 2 1] [2 2 2] [2 2 2] [2 2 3] [2 2 4]},...
    'fontsize',8,'gap',[0.07 0.07],'marg_w',[0.01 0.01],'marg_h',[0.04 0.04]);
axis([f1.s1 f1.s2 f1.s3 f1.s4 f1.s5 f1.s6],'image');

surf(f1.s1,test_I);% show intensity surface plot

plot(f1.s2,[0 size(test_I,2)],[0 size(test_I,1)],'linestyle','none');
plot(f1.s2,centroids(:,1),centroids(:,2),'rx');
plot(f1.s2,centerA(1),centerA(2),'r+','markersize',10);
plot(f1.s2,AL_grid2(:,1),AL_grid2(:,2),'go');%show the rotated grid
plot(f1.s2,size(test_I,2)/2,size(test_I,1)/2,'mo','markersize',16);

surf(f1.s3,test_I);% show intensity surface plot
plot(f1.s4,centroids(:,1),centroids(:,2),'rx','visible','off');
plot(f1.s4,centerA(1),centerA(2),'r+','markersize',10,'visible','off');
plot(f1.s4,AL_grid2(:,1),AL_grid2(:,2),'go','visible','off');
quiver(f1.s4,AL_grid2(:,1),AL_grid2(:,2),DX,DY,1,'color','y',...
    'maxheadsize',0.5);

scatter3(f1.s5,AL_grid2(:,1),AL_grid2(:,2),DX,'filled','cdata',DX);
axes(f1.s5);
f1.DX_h=plot(sfx);

scatter3(f1.s6,AL_grid2(:,1),AL_grid2(:,2),DY,'filled','cdata',DY);
axes(f1.s6);
f1.DY_h=plot(sfy);


% Plot and figure formatting
xylabels(f1.s1,'x','y');
xylabels(f1.s4,'x','y');
xylabels(f1.s5,'x','y');
xylabels(f1.s6,'x','y');
title(f1.s2,'Grid alignment','fontsize',12);
title(f1.s4,'Total distortion displacement field','fontsize',12);
title(f1.s5,'X distortion displacement field','fontsize',12);
title(f1.s6,'Y distortion displacement field','fontsize',12);
colormap(f1.s1,'gray');
colormap(f1.s3,'gray');
set(f1.s2,'ylim',f1.s1.YLim,'xlim',f1.s1.XLim,'position',f1.s1.Position);
set(f1.s4,'ylim',f1.s3.YLim,'xlim',f1.s3.XLim,'position',f1.s3.Position);
set(f1.s5,'dataaspectratio',[1 1 0.1]);
set(f1.s6,'dataaspectratio',[1 1 0.1]);
set(findall(f1.f,'type','axes'),'ydir','reverse');
linkaxes([f1.s1 f1.s2 f1.s3 f1.s4 f1.s5 f1.s6],'xy');
f1.f.Name=I1.file;

% Store additional fields in f2 structure output variable
f1.I1=I1;
f1.optimal_rot=optimal_rot;
f1.DX=DX;
f1.DY=DY;
f1.AL_grid2=AL_grid2;
f1.sfx=sfx;
f1.sfy=sfy;

disp('Distortion and grid alignment analysis completed');