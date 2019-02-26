% This fcn uses a mask to remove noise in the image
function [Iout]=mask_denoise(Iin)
%creak mask to remove noise
    I_mask=medfilt2(Iin,[3 3]);
    
    %binarize mask
    I_mask(~isnan(I_mask))=1;
    I_mask(isnan(I_mask))=0;
    
    % morphological operations
    I_mask=bwmorph(I_mask,'open');
    
    I_mask=abs(I_mask-1);%invert binary image
    I_mask=bwareaopen(I_mask,100);%remove small specks
    
    I_mask=bwmorph(I_mask,'open');
    I_mask=abs(I_mask-1);%invert the inverted image
    I_mask=double(I_mask);%convert to double array
    I_mask(I_mask==0)=nan;
    Iout=I_mask.*Iin;%remove noise