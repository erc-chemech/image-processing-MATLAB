function I_corr=rgb_correction(I,ref,varargin)
% Author: Joshua Yeh
% Date created: 2017/11/01
% 
%% DESCRIPTION
% This function performs an RGB correction to an rgb image, I, in reference
% to a white background. The input images should have values between 0 and
% 255.
% NOTE on 'simple' algorithm: Since the correction is based on a scaling
% factor determined by the mean of the pixels values of the white standard
% image, there will be a distribution of pixels values centerd at 255, when
% the image is naively corrected. However, in this case, pixel values
% greater tahn 255 is unphysical (in terms of grayscaling). Thus, pixels
% are corrected to a threshold value (~250) so that that there will not be
% an abundance of oversatuated pixels. Using this adjusment, any pixel
% values above 255 is replaced with values of 255 (saturation).
% 
%% SYNTAX
% I_corr=rgb_correction(I,ref,type,thresh,flag)
% 
%% INPUT VARIABLES
% 
% I: a 2D RGB image that will be color-corrected to a reference white
% background (double)
% ref: a reference (sub) image representing the white reference standard
% (double)
% type (optional): defines the type of rgb correction to be performed
%   -'simple': A simple color correction by using a reference white area
%   and assuming that r=b=g
% thresh (optional): threshold values in which the pixels values will be corrected to
% (default is 250)
% flag (optional): a flag that tells the program to replace values >255 with 255
% (saturation) when flag=1 (default flag is 1)
%
%% OUTPUT VARIABLES
% I_corr: corrected rgb image on normalized rgb ratio (Grassman's Law)
%
%%
% Determine number of inputs
switch nargin
    case 2
        type='simple';
        thresh=250;
        flag=1;        
    case 3
        type=varargin{1};
        thresh=250;
        flag=1;        
    case 4
        type=varargin{1};
        thresh=varargin{2};
        flag=1;
    case 5
        type=varargin{1};
        thresh=varargin{2};
        flag=varargin{3};
end

% convert I and ref to double
I=double(I);
ref=double(ref);
switch type
    case 'simple'
        %calibrate pixels to a white region by rescaling the pixels to
        %thresh
        for dum=1:3
            white=ref(:,:,dum);
            channel=I(:,:,dum);
            pd=fitdist(white(:),'normal');%fit normal dist to histogram
            
            % Determine offset factor
            offset=thresh./pd.mu;
            
            % Apply correction
            I_corr(:,:,dum)=channel.*offset;
            if flag==1
                ii=find(I_corr(:)>255);
                I_corr(ii)=255;
                disp(['# of saturated pixels detected: ',num2str(length(ii))]);
            end
        end
end
