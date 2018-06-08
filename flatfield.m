function out=flatfield(filenames,varargin)
% Author: Joshua Yeh
% 18/06/07
%% DESCRIPTION
% This script outputs a surface fitted intensity image of calibration
% images obtained from the Nikon AZ100 microscope.
% 
%% INPUT
% filenames: list of filenames
% 
% varargout: 'fieldname', <key>
    % 'in_size': pixel dimensions of input images
    % 'med_span': median filter square span size
% 
%% OUTPUT
% out: the output surface fitted intensity image
% 
%% PARSE THE INPUTS
narginchk(1,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('in_size',[512 512],@(x) isnumeric(x));
params.addParameter('med_span',5,@(x) isnumeric(x)&x>0&x<1);
params.parse(varargin{:});

% Extract out parameters from parsed input
in_size=params.Results.in_size;
med_span=params.Results.med_span;

%% Import images from filenames and perform statistical calculations

Isum=zeros(in_size);
for dum=1:numel(filenames)
    filename=filenames{dum};%get current filename
    [~,name,~]=fileparts(filename);%parse filename
    I.(name)=import_tiff_stack(filename,1);%Import tif file
    Isum=Isum+I.(name).tiff_stack;%add all of the images
end

% Compute the mean of the intensites
Imean=Isum./numel(filenames);
I_n=Imean./max(Imean(:));%normalize the intensities

% nominal image resolution
res=I.(name).info.UnknownTags(2).Value;%res in um/px
[x,y]=meshgrid(1:in_size(2),1:in_size(1));
a=linspace(min(x(:)),max(x(:)),in_size(1));
w=1;

[I_n2,A,B]=coord2image(x(:),y(:),I_n(:),w,'mean');
out=medfilt2(I_n2,[med_span med_span]);

%% PLOT THE RESULTS

f1=my_fig(1);
axis(f1.s1,'image');
set(f1.s1,'ydir','reverse');

imagesc(f1.s1,out);
xylabels(f1.s1,'x (px)','y (px)');
center_axes(f1.s1,'margins',10);


