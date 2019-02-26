function stats=DA_ROI(I1,X,Y,ROIs,varargin)
% This function performs image quantification of corrected confocal images
% of samples containing DA mechanophore (bond breakage). This function can
% handle multiple ROIs within the same image.

%% Parse input variables
narginchk(3,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('combine_px',0,@(x) x==0|x==1|islogical(x));
params.parse(varargin{:});

c_px=params.Results.combine_px;


for dum0=1:numel(ROIs)
    name=(['stats',num2str(dum0)]);%field name based on ROI index

    %extract subimage
    [ROI,subI,subx,suby]=ex_subI(I1,X,Y,ROIs{dum0});

    %binarize subimage
    bw=subI;    bw(~isnan(bw))=true;    bw(isnan(bw))=false;

    % perform morphological changes to the subimage
    bw=bwmorph(bw,'bridge');
    bw=bwmorph(bw,'fill');
    bw=bwmorph(bw,'open');

    %extract edge points
    bw=edge(bw);
    stats.(name)=regionprops(bw,'pixellist');
    
    if c_px==1% check to see if toggle to combine edge px lists is true | 1
        pl=[];
        for ddum=1:numel(stats.(name))
            cpl=stats.(name)(ddum).PixelList;%"current" pixel list
            if size(cpl,1)>100%ignore mislabeled edges
                pl=cat(1,pl,cpl);
            end
        end
        stats.(name)(1).PixelList=pl;
    end

    %find indices assoc. with being an actual number
    idx0=find(~isnan(subI(:)));
    [r0,c0]=ind2sub(size(subI),idx0);%find coord. assoc. w/ idx0
    subI_coord=subI(idx0);%intensities assoc. w/ idx0

    % find min distance of each pixel to the contour length
    [idx1,dist]=knnsearch(...
        stats.(name)(1).PixelList,[c0,r0]);
    stats.(name)(1).u_dist=unique(dist);

    % calculate the stats based on dist from edge in pxs
    stats.(name)(1).cummed_I=nan(numel(stats.(name)(1).u_dist),1);
    stats.(name)(1).cummean_I=nan(numel(stats.(name)(1).u_dist),1);
    stats.(name)(1).cumstd_I=nan(numel(stats.(name)(1).u_dist),1);
    stats.(name)(1).med_I=nan(numel(stats.(name)(1).u_dist),1);
    stats.(name)(1).mean_I=nan(numel(stats.(name)(1).u_dist),1);
    stats.(name)(1).std_I=nan(numel(stats.(name)(1).u_dist),1);
    for dum1=1:numel(stats.(name)(1).u_dist)
        %cumulative median
        stats.(name)(1).cummed_I(dum1)=...
            nanmedian(subI_coord(dist<=stats.(name)(1).u_dist(dum1)));

        %cumulative mean
        stats.(name)(1).cummean_I(dum1)=...
            nanmean(subI_coord(dist<=stats.(name)(1).u_dist(dum1)));

        %cumulative std
        stats.(name)(1).cumstd_I(dum1)=...
            nanstd(subI_coord(dist<=stats.(name)(1).u_dist(dum1)));

        %median
        stats.(name)(1).med_I(dum1)=...
            nanmedian(subI_coord(dist==stats.(name)(1).u_dist(dum1)));

        %mean
        stats.(name)(1).mean_I(dum1)=...
            nanmean(subI_coord(dist==stats.(name)(1).u_dist(dum1)));

        %std
        stats.(name)(1).std_I(dum1)=...
            std(subI_coord(dist==stats.(name)(1).u_dist(dum1)));

    end
end