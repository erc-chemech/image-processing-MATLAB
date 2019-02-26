% This fcn handles flat field correction and distortion corrections
function [Iout,X,Y]=corrections(Iin,IH,dx,dy,meta)
    offset=meta.offset;
    T1=meta.T1;
    w2=meta.w2;
    res=meta.res;

    %flat field correction
    plane_initial=(Iin-offset)./IH;%flat field correction
    plane_initial(IH(:)<T1)=0;%remove overcorrect areas, defined by T1
    
    % Convert px coord to cartesian coord
    [x,y]=meshgrid(1:size(plane_initial,2),1:size(plane_initial,1));
    x=x.*res;    y=y.*res; %convert to real position in um
    z1=plane_initial(:);%intensity values
    
    % Perform distortion correction
    x=x-dx;    y=y-dy;
    
    % Turn 2d array into column array
    [x,y]=prepareCurveData(x,y);
    
    % Performing binning of current image
    [Iout,X,Y]=coord2image(x,y,z1,w2,'mean');
    
end

