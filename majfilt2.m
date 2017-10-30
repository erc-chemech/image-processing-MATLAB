function I_out=majfilt2(varargin)
% Author: Joshua Yeh
% Date created: 2017/10/30
%
%% DESCRIPTION
% This functions applies a MAJORITY filter on an image using a user-defined
% square window.
% 
%% INPUT VARIABLES
% I: Image to be filtered (must be a 2d binary image, logical array)
% pad: scalar value that is used to define the number of pixels to the
%   left, right, top, and bottom of the center pixel that will define the
%   window in which the MAJORITY filter will be applied to. For example, if
%   pad=2, this will result in a square window of 5x5:
%   O O O O O
%   O O O O O
%   O O X O O
%   O O O O O
%   O O O O O
%   A 5x5 window with pad=2, where "X" denotes the center pixel.
%   The border of the image will be padded with nans. This border is defined by
%   pad.
% n: number of times to repeat the filter
%% OUTPUT VARIABLES
% I_out: Filtered image
% 
%%
% Determine number of input arguments
switch nargin
    case 2
        I=varargin{1};
        pad=varargin{2};
        n=1;
    case 3
        I=varargin{1};
        pad=varargin{2};
        n=varargin{3};
end

% Get dimensons on the image, I
[r,c]=size(I);

% Convert logical array to double
I=double(I);

% Copy image, I, without border defined by pad
I_inner=I;
I2=I;

% define range in which filter is applied
r_range=[pad+1,r-pad];%row range
c_range=[pad+1,c-pad];%column range

% Apply filter recursively on defined window (defined by pad)
for dum=1:n%apply filter n times
    for center_r=r_range(1):r_range(2)%iterate over rows
        for center_c=c_range(1):c_range(2)%iterate over columns
            sub_I=I2(center_r-pad:center_r+pad,center_c-pad:center_c+pad);
            I_inner(center_r,center_c)=mode(sub_I(:));
        end
    end
    I2=I_inner;
end

% Replace border of I_out with zeros (border thicknes is defined by pad)
I_out=nan(r,c);
I_out(pad+1:r-pad,pad+1:c-pad)=I_inner(pad+1:r-pad,pad+1:c-pad);
% I_out=logical(I_out);%convert array back to logical
