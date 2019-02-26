% This fcn extracts a subimage from the original image
function [ROI,subI,subx,suby]=ex_subI(Iin,X,Y,ROI)
    subI=Iin(ROI(1):ROI(2),ROI(3):ROI(4));%reg. of interest of the subimage
    subx=X(ROI(1):ROI(2),ROI(3):ROI(4));
    suby=Y(ROI(1):ROI(2),ROI(3):ROI(4));
