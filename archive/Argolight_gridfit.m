function f2=Argolight_gridfit(I1,test,varargin)
% Author: Joshua Yeh
% Date created: 2018/01/26
%% DESCRIPTION
% This script is only for tiff images taken with Nikon AZ100! Also, it is
% important that the calibration image is not severely rotated (no more
% than +/- 10 degrees in rotation).
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
% test: the image array in which the fir will be based on
% 
% T: threshold value for binarizing image (for identification of grid
% points)
% 
% type: grid fitting algorithm
    % 'cross': coarse fitting based on the cross mark
    % 'points': fine-tuned fitting based on minimizing the total displacement
    % magnitudes of all the grid points (this algorithm is based on the
    % cross algorithm)
%%

% Parse input variables
narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('T',100,@(x) isnumeric(x));
% params.addParameter('type','points',@(x) ischar(x));
params.parse(varargin{:});

% determine what type of rotation algorithm to use
% if strcmp(params.type,'cross')
%     flag=1;
% elseif strcmp(params.type,'points')
%     flag=2;
% end

% create anonymous fcns
% xy rotation matrix (counterclockwise) in degs
rotxy=@(x) [cosd(x) -sind(x);sind(x) cosd(x)];

% Binarize the average image
I1.bw=imbinarize(test,params.Results.T);
I1.bw=bwmorph(I1.bw,'fill');

% Extract the pixel resolution (um/px)
px_res=I1.info(1).UnknownTags(2).Value;
disp(['res: ',num2str(px_res),' ',956,'m/px']);

% Find location of grid points
I1.stats=regionprops(I1.bw,'centroid','area','image','PixelList');
centroids=reshape([I1.stats.Centroid],[2,size(I1.stats,1)])';

% % Find the cross point (if any) by identifying outliers in the grid point
% % area and the centroids of identified grid points
% TF=isoutlier([I1.stats(:).Area]);
% iio=find(TF==1);
% [~,iioa]=max([I1.stats(iio).Area]);
% iio=iio(iioa);
% cross=[I1.stats(iio).Centroid];
% cross_cr=I1.stats(iio).PixelList;
% cross_subimage=test(min(cross_cr(:,2)):...
%     max(cross_cr(:,2)),...
%     min(cross_cr(:,1)):...
%     max(cross_cr(:,1)));


% Identify an approximate middle grid point
midx=size(test,2)/2;
midy=size(test,1)/2;
D=sqrt((centroids(:,1)-midx).^2+(centroids(:,2)-midy).^2);
[~,m3]=min(D);
centerA=centroids(m3,:);

% Construct grid assuming no rotation
spacing=50/px_res;%spacing between dots in pixels
R_row=centerA(1):spacing:I1.info(1).Width;
L_row=centerA(1):-spacing:1;
row=sort([L_row R_row],'ascend');
L_col=centerA(2):spacing:I1.info(1).Height;
U_col=centerA(2):-spacing:1;
col=sort([U_col L_col],'ascend');

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

% % Perform coarse rotation fit based on cross marker
% IX=I1.stats(iio).Image;
% 
% % First create model cross
% TX=zeros(size(IX));
% 
% % add horizontal striping
% if mod(size(TX,1),2)==0
%     TX(size(TX,1)/2:size(TX,1)/2+1,3:end-2)=ones(2,size(TX,2)-4);
% else
%     TX(ceil(size(TX,1)/2),:)=ones(1,size(TX,2));
% end
% 
% % add vertical striping
% if mod(size(TX,2),2)==0
%     TX(3:end-2,size(TX,2)/2:size(TX,2)/2+1)=ones(size(TX,1)-4,2);
% else
%     TX(3:end-2,ceil(size(TX,2)/2))=ones(size(TX,1)-4,1);
% end
% 
% % Rotate model cross and calculate product of two images
% rots=-45:0.05:45;
% eX=nan(length(rots),1);
% count=1;
% mX=0;
% 
% % cycle through a range of rotations defined by rots (rotations)
% for rot=rots
%     rot_TX=imrotate(TX,rot);
% 
%     %Check image size match, otherwise pad images so that size matches
%     if sum(size(rot_TX)==size(IX))==2
%         mX=rot_TX.*IX;%multiple images
%         mX=sum(mX(:));
% 
%     else
%         padsize1=(size(rot_TX,1)-size(IX,1))/2;%padding for top, bottom
%         padsize2=(size(rot_TX,2)-size(IX,2))/2;%padding for left, right
% 
%         if mod(padsize1,1)==1&&mod(padsize2,1)==1
%             pad_IX=[zeros(padsize1,size(IX,2));IX;zeros(padsize1,size(IX,2))];
%             pad_IX=[zeros(size(pad_IX,1),padsize2),pad_IX,...
%                 zeros(size(pad_IX,1),padsize2)];
%         elseif mod(padsize1,1)~=0&&mod(padsize2,0)==1&&padsize1>1
%             % more top padding
%             pad_IX_top=[zeros(padsize1+0.5,size(IX,2));IX;...
%                 zeros(padsize1-0.5,size(IX,2))];
%             pad_IX_top=[zeros(size(pad_IX_top,1),padsize2),pad_IX_top,...
%                 zeros(size(pad_IX_top,1),padsize2)];
% 
%             % more bottom padding
%             pad_IX_bottom=[zeros(padsize1-0.5,size(IX,2));IX;...
%                 zeros(padsize1+0.5,size(IX,2))];
%             pad_IX_bottom=[zeros(size(pad_IX_bottom,1),padsize2),...
%                 pad_IX_bottom,zeros(size(pad_IX_bottom,1),padsize2)];
% 
%             % find best fit
%             mX_top=rot_TX.*pad_IX_top;
%             mX_top=sum(mX_top(:));
%             mX_bottom=rot_TX.*pad_IX_bottom;
%             mX_bottom=sum(mX_bottom(:));
%             mX=max([mX_top mX_bottom]);
% 
%         elseif mod(padsize1,1)~=0&&mod(padsize2,1)==0&&padsize1<1
%             % padding ontop
%             pad_IX_top=[zeros(1,size(IX,2));IX];
%             pad_IX_top=[zeros(size(pad_IX_top,1),padsize2),pad_IX_top,...
%                 zeros(size(pad_IX_top,1),padsize2)];
% 
%             % padding on the bottom
%             pad_IX_bottom=[IX;zeros(1,size(IX,2))];
%             pad_IX_bottom=[zeros(size(pad_IX_bottom,1),padsize2),...
%                 pad_IX_bottom,zeros(size(pad_IX_bottom,1),padsize2)];
% 
%             % find best fit
%             mX_top=rot_TX.*pad_IX_top;
%             mX_top=sum(mX_top(:));
%             mX_bottom=rot_TX.*pad_IX_bottom;
%             mX_bottom=sum(mX_bottom(:));
%             mX=max([mX_top mX_bottom]);
% 
%         elseif mod(padsize1,1)==0&&mod(padsize2,1)~=0&&padsize2>1
%             % more padding on the left
%             pad_IX_left=[zeros(size(IX,1),padsize2+0.5),IX,...
%                 zeros(size(IX,1),padsize2-0.5)];
%             pad_IX_left=[zeros(padsize1,size(pad_IX_left,2));pad_IX_left;...
%                 zeros(padsize1,size(pad_IX_left,2))];
% 
%             % more padding on the right
%             pad_IX_right=[zeros(size(IX,1),padsize2-0.5),IX,...
%                 zeros(size(IX,1),padsize2+0.5)];
%             pad_IX_right=[zeros(padsize1,size(pad_IX_right,2));pad_IX_right;...
%                 zeros(padsize1,size(pad_IX_right,2))];
% 
%             % find best fit
%             mX_left=rot_TX.*pad_IX_left;
%             mX_left=sum(mX_left(:));
%             mX_right=rot_TX.*pad_IX_right;
%             mX_right=sum(mX_right(:));
%             mX=max([mX_left mX_right]);
% 
%         elseif mod(padsize1,1)==0&&mod(padsize2,1)~=0&&padsize2<1
%             % more padding on the left
%             pad_IX_left=[zeros(size(IX,1),1),IX];
%             pad_IX_left=[zeros(padsize1,size(pad_IX_left,2));pad_IX_left;...
%                 zeros(padsize1,size(pad_IX_left,2))];
% 
%             % more padding on the right
%             pad_IX_right=[IX,zeros(size(IX,1),1)];
%             pad_IX_right=[zeros(padsize1,size(pad_IX_right,2));pad_IX_right;...
%                 zeros(padsize1,size(pad_IX_right,2))];
% 
%             % find best fit
%             mX_left=rot_TX.*pad_IX_left;
%             mX_left=sum(mX_left(:));
%             mX_right=rot_TX.*pad_IX_right;
%             mX_right=sum(mX_right(:));
%             mX=max([mX_left mX_right]);
% 
%         elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
%                 padsize1>1&&padsize2>1
%             % more padding on the left and top
%             pad_IX_LT=[zeros(size(IX,1),padsize2+0.5),IX,...
%                 zeros(size(IX,1),padsize2-0.5)];
%             pad_IX_LT=[zeros(padsize1+0.5,size(pad_IX_LT,2));pad_IX_LT;...
%                 zeros(padsize1-0.5,size(pad_IX_LT,2))];
% 
%             % more padding on the left and bottom
%             pad_IX_LB=[zeros(size(IX,1),padsize2+0.5),IX,...
%                 zeros(size(IX,1),padsize2-0.5)];
%             pad_IX_LB=[zeros(padsize1-0.5,size(pad_IX_LB,2));pad_IX_LB;...
%                 zeros(padsize1+0.5,size(pad_IX_LB,2))];
% 
%             % more padding on the right and top
%             pad_IX_RT=[zeros(size(IX,1),padsize2-0.5),IX,...
%                 zeros(size(IX,1),padsize2+0.5)];
%             pad_IX_RT=[zeros(padsize1+0.5,size(pad_IX_RT,2));pad_IX_RT;...
%                 zeros(padsize1-0.5,size(pad_IX_RT,2))];
% 
%             % more padding on the right and bottom
%             pad_IX_RB=[zeros(size(IX,1),padsize2-0.5),IX,...
%                 zeros(size(IX,1),padsize2+0.5)];
%             pad_IX_RB=[zeros(padsize1-0.5,size(pad_IX_RB,2));pad_IX_RB;...
%                 zeros(padsize1+0.5,size(pad_IX_RB,2))];
% 
%             % find best fit
%             mX_IX_LT=rot_TX.*pad_IX_LT;
%             mX_IX_LT=sum(mX_IX_LT(:));
%             mX_IX_LB=rot_TX.*pad_IX_LB;
%             mX_IX_LB=sum(mX_IX_LB(:));
%             mX_IX_RT=rot_TX.*pad_IX_RT;
%             mX_IX_RT=sum(mX_IX_RT(:));
%             mX_IX_RB=rot_TX.*pad_IX_RB;
%             mX_IX_RB=sum(mX_IX_RB(:));
%             mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
% 
%         elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
%                 padsize2<1&&padsize1>1
%             % more padding on the left and top
%             pad_IX_LT=[zeros(size(IX,1),1),IX];
%             pad_IX_LT=[zeros(padsize1+0.5,size(pad_IX_LT,1));pad_IX_LT;...
%                 zeros(padsize1-0.5,size(pad_IX_LT,1))];
% 
%             % more padding on the left and bottom
%             pad_IX_LB=[zeros(size(IX,1),1),IX];
%             pad_IX_LB=[zeros(padsize1-0.5,size(pad_IX_LB,1));pad_IX_LB;...
%                 zeros(padsize1+0.5,size(pad_IX_LB,1))];
% 
%             % more padding on the right and top
%             pad_IX_RT=[IX,zeros(size(IX,1),1)];
%             pad_IX_RT=[zeros(padsize1+0.5,size(pad_IX_RT,1));pad_IX_RT;...
%                 zeros(padsize1-0.5,size(pad_IX_RT,1))];
% 
%             % more padding on the right and bottom
%             pad_IX_RB=[IX,zeros(size(IX,1),1)];
%             pad_IX_RB=[zeros(padsize1-0.5,size(pad_IX_RB,1));pad_IX_RB;...
%                 zeros(padsize1+0.5,size(pad_IX_RB,1))];
% 
%             % find best fit
%             mX_IX_LT=rot_TX.*pad_IX_LT;
%             mX_IX_LT=sum(mX_IX_LT(:));
%             mX_IX_LB=rot_TX.*pad_IX_LB;
%             mX_IX_LB=sum(mX_IX_LB(:));
%             mX_IX_RT=rot_TX.*pad_IX_RT;
%             mX_IX_RT=sum(mX_IX_RT(:));
%             mX_IX_RB=rot_TX.*pad_IX_RB;
%             mX_IX_RB=sum(mX_IX_RB(:));
%             mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
% 
%         elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
%                 padsize2>1&&padsize1<1
%             % more padding on the left and top
%             pad_IX_LT=[zeros(size(IX,1),padsize2+0.5),IX,...
%                 zeros(size(IX,1),padsize2-0.5)];
%             pad_IX_LT=[zeros(1,size(pad_IX_LT,2));pad_IX_LT];
% 
%             % more padding on the left and bottom
%             pad_IX_LB=[zeros(size(IX,1),padsize2+0.5),IX,...
%                 zeros(size(IX,1),padsize2-0.5)];
%             pad_IX_LB=[pad_IX_LB;zeros(1,size(pad_IX_LB,2))];
% 
%             % more padding on the right and top
%             pad_IX_RT=[zeros(size(IX,1),padsize2-0.5),IX,...
%                 zeros(size(IX,1),padsize2+0.5)];
%             pad_IX_RT=[zeros(1,size(pad_IX_RT,2));pad_IX_RT];
% 
%             % more padding on the right and bottom
%             pad_IX_RB=[zeros(size(IX,1),padsize2-0.5),IX,...
%                 zeros(size(IX,1),padsize2+0.5)];
%             pad_IX_RB=[pad_IX_RB;zeros(1,size(pad_IX_RB,2))];
% 
%             % find best fit
%             mX_IX_LT=rot_TX.*pad_IX_LT;
%             mX_IX_LT=sum(mX_IX_LT(:));
%             mX_IX_LB=rot_TX.*pad_IX_LB;
%             mX_IX_LB=sum(mX_IX_LB(:));
%             mX_IX_RT=rot_TX.*pad_IX_RT;
%             mX_IX_RT=sum(mX_IX_RT(:));
%             mX_IX_RB=rot_TX.*pad_IX_RB;
%             mX_IX_RB=sum(mX_IX_RB(:));
%             mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
% 
%         elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
%                 padsize1<1&&padsize2<1
%             % more padding on the left and top
%             pad_IX_LT=[zeros(size(IX,1),1),IX];
%             pad_IX_LT=[zeros(1,size(pad_IX_LT,2));pad_IX_LT];
% 
%             % more padding on the left and bottom
%             pad_IX_LB=[zeros(size(IX,1),1),IX];
%             pad_IX_LB=[pad_IX_LB;zeros(1,size(pad_IX_LB,2))];
% 
%             % more padding on the right and top
%             pad_IX_RT=[IX,zeros(size(IX,1),1)];
%             pad_IX_RT=[zeros(1,size(pad_IX_RT,2));pad_IX_RT];
% 
%             % more padding on the right and bottom
%             pad_IX_RB=[IX,zeros(size(IX,1),1)];
%             pad_IX_RB=[pad_IX_RB;zeros(1,size(pad_IX_RB,2))];
% 
%             % find best fit
%             mX_IX_LT=rot_TX.*pad_IX_LT;
%             mX_IX_LT=sum(mX_IX_LT(:));
%             mX_IX_LB=rot_TX.*pad_IX_LB;
%             mX_IX_LB=sum(mX_IX_LB(:));
%             mX_IX_RT=rot_TX.*pad_IX_RT;
%             mX_IX_RT=sum(mX_IX_RT(:));
%             mX_IX_RB=rot_TX.*pad_IX_RB;
%             mX_IX_RB=sum(mX_IX_RB(:));
%             mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
% 
%         end
%     end
%     eX(count)=mX;
%     count=count+1;%update counter
% end
% ee1=find(eX==max(eX));
% ee2=round(length(ee1)/2);
% optimal_rot=rots(ee1(ee2));
% disp(['Optimal rotation based of cross fit (deg.): ',num2str(rots(ee1(1))),...
%     ' to ',num2str(rots(ee1(end)))]);

% finer rotation fit based on coarse rotation fit based on cross fit
% if flag==2
    rots=linspace(-5,5,200);
    e_rots=[];
    for rot=rots
        % locate corresponding points
        pt_store=[];
        for dum=1:size(AL_grid,2)
            AL_grid_rot(dum,:)=(AL_grid(dum,:)-centerA)*rotxy(rot)+centerA;
            
            %find nearest point
            dx=centroids(:,1)-AL_grid_rot(dum,1);
            dy=centroids(:,2)-AL_grid_rot(dum,2);
            dtot=sqrt(dx.^2+dy.^2);
            [~,ii]=min(dtot);
            pt_store=[pt_store;ii dtot(ii)];
        end
        e_rots=[e_rots;nansum(pt_store(:,2))];%store total displacement
    end
    [~,ii]=min(e_rots);%rotation w/ least total displacement is optimal fit
    optimal_rot=rots(ii);
    disp(['Optimal rotation (deg.): ',num2str(optimal_rot)]);
% end


% Rotate grid by the determined optimal rotation and recenter based on
% cross location
for dum=1:size(AL_grid2,1)
    AL_grid2(dum,:)=AL_grid2(dum,:)*rotxy(optimal_rot)+centerA;
end

% Determine the x and y displacement of the reference rotated grid and the
% imaged grid patter
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


%% Fit a surface to x displacement field

disp('Performing x displacement surface fit');
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';

[sfx, gof] = fit( [AL_grid2(:,1), AL_grid2(:,2)], DX, ft, opts );

%% Fit a surface to y displacement field

disp('Performing y displacement surface fit');
[sfy, gof] = fit( [AL_grid2(:,1), AL_grid2(:,2)], DY, ft, opts );


%% PLOT RESULTS
disp('Plotting results...');
f2=my_fig(99,{[2 2 1] [2 2 1] [2 2 2] [2 2 2] [2 2 3] [2 2 4]},...
    'fontsize',8,'gap',[0.07 0.07],'marg_w',[0.01 0.01],'marg_h',[0.04 0.04]);
axis([f2.s1 f2.s2 f2.s3 f2.s4 f2.s5 f2.s6],'image');

surf(f2.s1,test);% show intensity surface plot

plot(f2.s2,[0 size(test,2)],[0 size(test,1)],'linestyle','none');
plot(f2.s2,centroids(:,1),centroids(:,2),'rx');
plot(f2.s2,centerA(1),centerA(2),'r+','markersize',10);
% plot(f2.s2,cross(1),cross(2),'g+','linewidth',2);
plot(f2.s2,AL_grid2(:,1),AL_grid2(:,2),'go');%show the rotated grid

surf(f2.s3,test);% show intensity surface plot
plot(f2.s4,centroids(:,1),centroids(:,2),'rx','visible','off');
plot(f2.s4,centerA(1),centerA(2),'r+','markersize',10,'visible','off');
plot(f2.s4,AL_grid2(:,1),AL_grid2(:,2),'go','visible','off');
quiver(f2.s4,AL_grid2(:,1),AL_grid2(:,2),DX,DY,1,'color','y',...
    'maxheadsize',0.5);


scatter3(f2.s5,AL_grid2(:,1),AL_grid2(:,2),DX,'filled','cdata',DX);
axes(f2.s5);
f2.DX_h=plot(sfx,[AL_grid2(:,1), AL_grid2(:,2)], DX);

scatter3(f2.s6,AL_grid2(:,1),AL_grid2(:,2),DY,'filled','cdata',DY);
axes(f2.s6);
f2.DY_h=plot(sfy,[AL_grid2(:,1), AL_grid2(:,2)], DY);


% Plot and figure formatting
xylabels(f2.s1,'x','y');
xylabels(f2.s4,'x','y');
xylabels(f2.s5,'x','y');
xylabels(f2.s6,'x','y');
title(f2.s2,'Grid alignment','fontsize',12);
title(f2.s4,'Total distortion displacement field','fontsize',12);
title(f2.s5,'X distortion displacement field','fontsize',12);
title(f2.s6,'Y distortion displacement field','fontsize',12);
colormap(f2.s1,'gray');
colormap(f2.s3,'gray');
set(f2.s2,'ylim',f2.s1.YLim,'xlim',f2.s1.XLim,'position',f2.s1.Position);
set(f2.s4,'ylim',f2.s3.YLim,'xlim',f2.s3.XLim,'position',f2.s3.Position);
set(f2.s5,'dataaspectratio',[1 1 0.1]);
set(f2.s6,'dataaspectratio',[1 1 0.1]);
set([f2.DX_h f2.DY_h],'markerfacecolor','none','markeredgecolor','k');
set(findall(f2.f,'type','axes'),'ydir','reverse');
linkaxes([f2.s1 f2.s2 f2.s3 f2.s4 f2.s5 f2.s6],'xy');
f2.f.Name=I1.file;

% Store additional fields in f2 structure output variable
f2.I1=I1;
f2.optimal_rot=optimal_rot;
f2.DX=DX;
f2.DY=DY;
f2.AL_grid2=AL_grid2;
f2.sfx=sfx;
f2.sfy=sfy;