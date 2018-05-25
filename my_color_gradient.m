function out=my_color_gradient(color1,color2,n)
% Author: Joshua Yeh
% Date: 18/04/23
% Description: Creates a discretize gradient between two colors. This
% function is intended to create customized colormap for plotting purposes.
%% INPUT VARIABLES
% color1: rgb color vector for the left or bottom end of the colorbar. 
% Elements must be between 0 and 1.
%
% color2: rgb color vector for the right or top end of the colorbar.
% Elements must be between 0 and 1.
%
% n: number of discretize steps between the two colors
%% OUTPUt variables
% out: the colormap gradient
%
%%

% process the color1 variable
if ischar(color1)
    color1=char2color(color1);%if it is a matlab char ('r','g','b',etc.)
else
%     color1=color1./max(color1);% normalize colors
end

% process the color2 variable
if ischar(color2)
    color2=char2color(color2);%if it is a matlab char ('r','g','b',etc.)
else
%     color2=color2./max(color2);% normalize colors
end

% red channel
r_space=linspace(color1(1),color2(1),n);

% green channel
g_space=linspace(color1(2),color2(2),n);

% blue channel
b_space=linspace(color1(3),color2(3),n);

out=flipud([r_space;g_space;b_space]');
end

function out=char2color(color)
    switch color
        case 'r'
            out=[1 0 0];
        case 'g'
            out=[0 1 0];
        case 'b'
            out=[0 0 1];
        case 'k'
            out=[0 0 0];
        case 'y'
            out=[1 1 0];
        case 'c'
            out=[0 1 1];
        case 'm'
            out=[1 0 1];
        case 'w'
            out=[1 1 1];
    end
end