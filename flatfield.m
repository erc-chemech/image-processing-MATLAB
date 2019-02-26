function [IH,x0,y0,I_n]=flatfield(filenames,varargin)
% Author: Joshua Yeh
% 18/06/07
%% DESCRIPTION
% This script outputs a surface fitted intensity image of calibration
% images obtained from the Nikon AZ100 microscope.
% 
%% INPUT
% filenames: list of filenames (cell variable)
% 
% varargout: 'fieldname', <key>
    % 'in_size': pixel dimensions of input images
    % 'fit_sampling': value from 0 to 1 used to down sample image pts for
    % fitting
% 
%% OUTPUT
% out: the output surface fitted intensity image
% 
%% PARSE THE INPUTS
narginchk(1,inf);%check number of inputs is correct
params=inputParser;
params.CaseSensitive=false;
params.addParameter('in_size',[512 512],@(x) isnumeric(x));
params.addParameter('fit_sampling',0.1,@(x) isnumeric(x)&x>0&x<1);
params.parse(varargin{:});

% Extract out parameters from parsed input
in_size=params.Results.in_size;
fit_sampling=params.Results.fit_sampling;

%% Import images from filenames and perform statistical calculations

Isum=zeros(in_size);
for dum=1:numel(filenames)
    filename=filenames{dum};%get current filename
    [~,name,~]=fileparts(filename);%parse filename
    if isvarname(name)==0
        name='IH_name';
    end
    I.(name)=import_tiff_stack(filename,1);%Import tif file
    Isum=Isum+I.(name).tiff_stack;%add all of the images
end

% Compute the mean of the intensites
Imean=Isum./numel(filenames);
I_n=Imean./max(Imean(:));%normalize the intensities
[x0,y0]=meshgrid(1:size(I_n,2),1:size(I_n,1));

I_n2=imresize(I_n,fit_sampling);
linx=linspace(1,size(I_n,2),size(I_n2,2));
liny=linspace(1,size(I_n,1),size(I_n2,1));
[x,y]=meshgrid(linx,liny);

% Perform surface fitting
[xData, yData, zData] = prepareSurfaceData( x(:), y(:), I_n2(:) );

% Set up fittype and options.
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Span = 0.10;

% Fit model to data.
disp('Fitting...(takes awhile)');
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
IH=fitresult(x0,y0);
disp('Fitting done');

%% PLOT THE RESULTS

f1=my_fig(1);
axis(f1.s1,'image');
set(f1.s1,'ydir','reverse');

imagesc(f1.s1,IH);
xylabels(f1.s1,'x (px)','y (px)');
center_axes(f1.s1,'margins',10);


