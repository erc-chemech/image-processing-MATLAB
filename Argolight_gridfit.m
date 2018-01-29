function Argolight_gridfit(I1,test)
% Author: Joshua Yeh
% Date created: 2018/01/26
%% DESCRIPTION
% This script attempts to fit an Argolight field of rings pattern
% (Argorlight LM) to a calibration image This script also assumes that the
% tiff image that was imported was exported by the Nikon NIS Viewer. This
% is important in order to obtain accurate Metadata information. Grid
% points are spaced 50 microns apart.
% 
%% INPUT VARIABLES
% I1: structure variable that is exported using the import_tiff_stacks
% function
%%

% Binarize the average image
I1.bw=imbinarize(I1.ave_tiff,100);
I1.bw=bwmorph(I1.bw,'fill');

% Extract the pixel resolution (um/px)
px_res=I1.info(1).UnknownTags(2).Value;

% Find location of grid points
I1.stats=regionprops(I1.bw,'centroid','area','image');

% Find the cross point (if any) by identifying outliers in the grid point
% area and the centroids of identified grid points
TF=isoutlier([I1.stats(:).Area]);
iio=find(TF==1);
I1.stats2=I1.stats;
I1.stats2(iio)=[];
cross=[I1.stats(iio).Centroid];
centroids=reshape([I1.stats2.Centroid],[2,size(I1.stats2,1)])';

% Identify an approximate middle grid point
midx=size(I1.ave_tiff,2)/2;
midy=size(I1.ave_tiff,1)/2;
D=sqrt((centroids(:,1)-midx).^2+(centroids(:,2)-midy).^2);
[~,m3]=min(D);
centerA=centroids(m3,:);

% Estimate rotation based on cross
IX=I1.stats(iio).Image;

% First create model cross
TX=zeros(size(IX));

% add horizontal striping
if mod(size(TX,1),2)==0
    TX(size(TX,1)/2:size(TX,1)/2+1,3:end-2)=ones(2,size(TX,2)-4);
else
    TX(ceil(size(TX,1)/2),:)=ones(1,size(TX,2));
end

% add vertical striping
if mod(size(TX,1),2)==0
    TX(3:end-2,size(TX,2)/2:size(TX,2)/2+1)=ones(size(TX,1)-4,2);
else
    TX(:,ceil(size(TX,2)/2))=ones(size(TX,1),2);
end


% Rotate model cross and calculate product of two images
rots=0:90;
eX=nan(length(rots),1);
count=1;
for rot=rots
    rot_TX=imrotate(TX,rot);
    
    %Check image size match
    if sum(size(rot_TX)==size(IX))==2
        mX=rot_TX.*IX;%multiple images
        mX=sum(mX(:));
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START HERE%%%%%%%%%%%%
        padsize1=(size(rot_TX,1)-size(IX,1))/2;%padding for top, bottom
        padsize2=(size(rot_TX,2)-size(IX,2))/2;%padding for left, right
        
        if mod(padsize1,1)==1&&mod(padsize2,1)==1
            pad_IX=[zeros(padsize1,size(IX,2));IX;zeros(padsize1,size(IX,2))];
            pad_IX=[zeros(size(pad_IX,1),padsize2),pad_IX,...
                zeros(size(pad_IX,1),padsize2)];
        elseif mod(padsize1,1)~=0&&mod(padsize2,0)==1&&padsize1>1
            % more top padding
            pad_IX_top=[zeros(padsize1+0.5,size(IX,2));IX;...
                zeros(padsize1-0.5,size(IX,2))];
            pad_IX_top=[zeros(size(pad_IX_top,1),padsize2),pad_IX_top,...
                zeros(size(pad_IX_top,1),padsize2)];
            
            % more bottom padding
            pad_IX_bottom=[zeros(padsize1-0.5,size(IX,2));IX;...
                zeros(padsize1+0.5,size(IX,2))];
            pad_IX_bottom=[zeros(size(pad_IX_bottom,1),padsize2),...
                pad_IX_bottom,zeros(size(pad_IX_bottom,1),padsize2)];
            
            % find best fit
            mX_top=rot_TX.*pad_IX_top;
            mX_top=sum(mX_top(:));
            mX_bottom=rot_TX.*pad_IX_bottom;
            mX_bottom=sum(mX_bottom(:));
            mX=max([mX_top mX_bottom]);
            
        elseif mod(padsize1,1)~=0&&mod(padsize2,1)==0&&padsize1<1
            % padding ontop
            pad_IX_top=[zeros(1,size(IX,2));IX];
            pad_IX_top=[zeros(size(pad_IX_top,1),padsize2),pad_IX_top,...
                zeros(size(pad+IX_top,1),padsize2)];
            
            % padding on the bottom
            pad_IX_bottom=[IX;zeros(1,size(IX,2))];
            pad_IX_bottom=[zeros(size(pad_IX_bottom,1),padsize2),...
                pad_IX_bottom,zeros(size(pad_IX_bottom,1),padsize2)];
            
            % find best fit
            mX_top=rot_TX.*pad_IX_top;
            mX_top=sum(mX_top(:));
            mX_bottom=rot_TX.*pad_IX_bottom;
            mX_bottom=sum(mX_bottom(:));
            mX=max([mX_top mX_bottom]);
            
        elseif mod(padsize1,1)==0&&mod(padsize2,1)~=0&&padsize2>1
            % more padding on the left
            pad_IX_left=[zeros(size(IX,1),padsize2+0.5),IX,...
                zeros(size(IX,1),padsize2-0.5)];
            pad_IX_left=[zeros(padsize1,size(pad_IX_left,2));pad_IX_left;...
                zeros(padsize1,size(pad_IX_left,2))];
            
            % more padding on the right
            pad_IX_right=[zeros(size(IX,1),padsize2-0.5),IX,...
                zeros(size(IX,1),padsize2+0.5)];
            pad_IX_right=[zeros(padsize1,size(pad_IX_right,2));pad_IX_right;...
                zeros(padsize1,size(pad_IX_right,2))];
            
            % find best fit
            mX_left=rot_TX.*pad_IX_left;
            mX_left=sum(mX_left(:));
            mX_right=rot_TX.*pad_IX_right;
            mX_right=sum(mX_right(:));
            mX=max([mX_left mX_right]);
            
        elseif mod(padsize1,1)==0&&mod(padsize2,1)~=0&&padsize2<1
            % more padding on the left
            pad_IX_left=[zeros(size(IX,1),1),IX];
            pad_IX_left=[zeros(padsize1,size(pad_IX_left,2));pad_IX_left;...
                zeros(padsize1,size(pad_IX_left,2))];
            
            % more padding on the right
            pad_IX_right=[IX,zeros(size(IX,1),1)];
            pad_IX_right=[zeros(padsize1,size(pad_IX_right,2));pad_IX_right;...
                zeros(padsize1,size(pad_IX_right,2))];
            
            % find best fit
            mX_left=rot_TX.*pad_IX_left;
            mX_left=sum(mX_left(:));
            mX_right=rot_TX.*pad_IX_right;
            mX_right=sum(mX_right(:));
            mX=max([mX_left mX_right]);
            
        elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
                padsize1>1&&padsize2>1
            % more padding on the left and top
            pad_IX_LT=[zeros(size(IX,1),padsize2+0.5),IX,...
                zeros(size(IX,1),padsize2-0.5)];
            pad_IX_LT=[zeros(padsize2+0.5,size(pad_IX_LT,2));pad_IX_LT;...
                zeros(padsize2-0.5,size(pad_IX_LT,2))];
            
            % more padding on the left and bottom
            pad_IX_LB=[zeros(size(IX,1),padsize2+0.5),IX,...
                zeros(size(IX,1),padsize2-0.5)];
            pad_IX_LB=[zeros(padsize2-0.5,size(pad_IX_LB,2));pad_IX_LB;...
                zeros(padsize2+0.5,size(pad_IX_LB,2))];
            
            % more padding on the right and top
            pad_IX_RT=[zeros(size(IX,1),padsize2-0.5),IX,...
                zeros(size(IX,1),padsize2+0.5)];
            pad_IX_RT=[zeros(padsize2+0.5,size(pad_IX_RT,2));pad_IX_RT;...
                zeros(padsize2-0.5,size(pad_IX_RT,2))];
            
            % more padding on the right and bottom
            pad_IX_RB=[zeros(size(IX,1),padsize2-0.5),IX,...
                zeros(size(IX,1),padsize2+0.5)];
            pad_IX_RB=[zeros(padsize2-0.5,size(pad_IX_RB,2));pad_IX_RB;...
                zeros(padsize2+0.5,size(pad_IX_RB,2))];
            
            % find best fit
            mX_IX_LT=rot_TX.*pad_IX_LT;
            mX_IX_LT=sum(mX_IX_LT(:));
            mX_IX_LB=rot_TX.*pad_IX_LB;
            mX_IX_LB=sum(mX_IX_LB(:));
            mX_IX_RT=rot_TX.*pad_IX_RT;
            mX_IX_RT=sum(mX_IX_RT(:));
            mX_IX_RB=rot_TX.*pad_IX_RB;
            mX_IX_RB=sum(mX_IX_RB(:));
            mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
            
        elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
                padsize1<1&&padsize2>1
            % more padding on the left and top
            pad_IX_LT=[zeros(size(IX,1),1),IX];
            pad_IX_LT=[zeros(size(pad_IX_LT,1),padsize2+0.5);pad_IX_LT;...
                zeros(size(pad_IX_LT,1),padsize2-0.5)];
            
            % more padding on the left and bottom
            pad_IX_LB=[zeros(size(IX,1),1),IX];
            pad_IX_LB=[zeros(size(pad_IX_LB,1),padsize2-0.5);pad_IX_LB;...
                zeros(size(pad_IX_LB,1),padsize2+0.5)];
            
            % more padding on the right and top
            pad_IX_RT=[IX,zeros(size(IX,1),1)];
            pad_IX_RT=[zeros(size(pad_IX_RT,1),padsize2+0.5);pad_IX_RT;...
                zeros(size(pad_IX_RT,1),padsize2-0.5)];
            
            % more padding on the right and bottom
            pad_IX_RB=[IX,zeros(size(IX,1),1)];
            pad_IX_RB=[zeros(size(pad_IX_RB,1),padsize2-0.5);pad_IX_RB;...
                zeros(size(pad_IX_RB,1),padsize2+0.5)];
            
            % find best fit
            mX_IX_LT=rot_TX.*pad_IX_LT;
            mX_IX_LT=sum(mX_IX_LT(:));
            mX_IX_LB=rot_TX.*pad_IX_LB;
            mX_IX_LB=sum(mX_IX_LB(:));
            mX_IX_RT=rot_TX.*pad_IX_RT;
            mX_IX_RT=sum(mX_IX_RT(:));
            mX_IX_RB=rot_TX.*pad_IX_RB;
            mX_IX_RB=sum(mX_IX_RB(:));
            mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
            
        elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
                padsize1>1&&padsize2<1
            % more padding on the left and top
            pad_IX_LT=[zeros(size(IX,1),padsize2+0.5),IX,...
                zeros(size(IX,1),padsize2-0.5)];
            pad_IX_LT=[zeros(1,size(pad_IX_LT,2));pad_IX_LT];
            
            % more padding on the left and bottom
            pad_IX_LB=[zeros(size(IX,1),padsize2+0.5),IX,...
                zeros(size(IX,1),padsize2-0.5)];
            pad_IX_LB=[pad_IX_LB;zeros(1,size(pad_IX_LB,2))];
            
            % more padding on the right and top
            pad_IX_RT=[zeros(size(IX,1),padsize2-0.5),IX,...
                zeros(size(IX,1),padsize2+0.5)];
            pad_IX_RT=[zeros(1,size(pad_IX_RT,2));pad_IX_RT];
            
            % more padding on the right and bottom
            pad_IX_RB=[zeros(size(IX,1),padsize2-0.5),IX,...
                zeros(size(IX,1),padsize2+0.5)];
            pad_IX_RB=[pad_IX_RB;zeros(1,size(pad_IX_RB,2))];
            
            % find best fit
            mX_IX_LT=rot_TX.*pad_IX_LT;
            mX_IX_LT=sum(mX_IX_LT(:));
            mX_IX_LB=rot_TX.*pad_IX_LB;
            mX_IX_LB=sum(mX_IX_LB(:));
            mX_IX_RT=rot_TX.*pad_IX_RT;
            mX_IX_RT=sum(mX_IX_RT(:));
            mX_IX_RB=rot_TX.*pad_IX_RB;
            mX_IX_RB=sum(mX_IX_RB(:));
            mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
            
        elseif mod(padsize1,1)~=0&&mod(padsize2,1)~=0&&...
                padsize1<1&&padsize2<1
            % more padding on the left and top
            pad_IX_LT=[zeros(size(IX,1),1),IX];
            pad_IX_LT=[zeros(1,size(pad_IX_LT,2));pad_IX_LT];
            
            % more padding on the left and bottom
            pad_IX_LB=[zeros(size(IX,1),1),IX];
            pad_IX_LB=[pad_IX_LB;zeros(1,size(pad_IX_LB,2))];
            
            % more padding on the right and top
            pad_IX_RT=[IX,zeros(size(IX,1),1)];
            pad_IX_RT=[zeros(1,size(pad_IX_RT,2));pad_IX_RT];
            
            % more padding on the right and bottom
            pad_IX_RB=[IX,zeros(size(IX,1),1)];
            pad_IX_RB=[pad_IX_RB;zeros(1,size(pad_IX_RB,2))];
            
            % find best fit
            mX_IX_LT=rot_TX.*pad_IX_LT;
            mX_IX_LT=sum(mX_IX_LT(:));
            mX_IX_LB=rot_TX.*pad_IX_LB;
            mX_IX_LB=sum(mX_IX_LB(:));
            mX_IX_RT=rot_TX.*pad_IX_RT;
            mX_IX_RT=sum(mX_IX_RT(:));
            mX_IX_RB=rot_TX.*pad_IX_RB;
            mX_IX_RB=sum(mX_IX_RB(:));
            mX=max([mX_IX_LT mX_IX_LB mX_IX_RT mX_IX_RB]);
            
        end
    end
    eX(count)=mX;
    count=count+1;%update counter
end
ee1=find(eX==max(eX));
ee2=round(length(ee1)/2);
optimal_rot=-rots(ee1(ee2));

%% Construct grid assuming no rotation
spacing=50/px_res;%spacing between dots in pixels
R_row=cross(1):spacing:I1.info(1).Width;
L_row=cross(1):-spacing:1;
row=sort([L_row R_row],'ascend');
L_col=cross(2):spacing:I1.info(1).Height;
U_col=cross(2):-spacing:1;
col=sort([U_col L_col],'ascend');

AL_grid=nan(length(row)*length(col),2);
AL_grid2=nan(length(row)*length(col),2);
% AL_grid_I=zeros(length(row),length(col));
count=1;
for dum1=1:length(row)
    for dum2=1:length(col)
        AL_grid(count,:)=[row(dum1) col(dum2)];
        AL_grid2(count,:)=[row(dum1)-cross(1) col(dum2)-cross(2)];
        count=count+1;
    end
end

% Rotate grid by the determined optimal rotation
for dum=1:size(AL_grid2,1)
    AL_grid2(dum,:)=AL_grid2(dum,:)*[cosd(optimal_rot) -sind(optimal_rot);...
        sind(optimal_rot) cosd(optimal_rot)];
    AL_grid2(dum,:)=AL_grid2(dum,:)+cross;
end

f2=my_fig(99,{[1 1 1]},'fontsize',8);
axis(f2.s1,'image');

imagesc(f2.s1,test);
plot(f2.s1,centroids(:,1),centroids(:,2),'rx');
plot(f2.s1,centerA(1),centerA(2),'r+','markersize',10);
plot(f2.s1,cross(1),cross(2),'g+','linewidth',2);
plot(f2.s1,AL_grid(:,1),AL_grid(:,2),'gx');
plot(f2.s1,AL_grid2(:,1),AL_grid2(:,2),'go');

xylabels(f2.s1,'x','y');