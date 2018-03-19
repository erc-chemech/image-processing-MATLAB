function [x,y,varargout]=bin_data_array(X,Y,w,varargin)
%% DESCRIPTION
% This function sccepts a one-dimensional array and bins the data into bin
% size defined by a user-defined increment. The mean of the data inside
% each bin is calculated and outputted.
% 
%% INPUT VARIABLES
% X: the x column array
% Y: the y column array, must have the same size as X
% w: bin size in units of the X array variable
% 
% varargin: '<fieldname>',Value
    % 'ld': lower error bar (default is empty), if used, it must have the
    % same size as X
    % 
    % 'ud': upper error bar (default is empty), if used, it must have the
    % same size as Y
    %
    % 'binscale': bin scaling, user can choose between 'linear' and log'
    % (default is 'linear'). Logarithmic scaling is base 10.
% 
%% OUTPUT VARIABLES
% x: bin y array
% y: binned y array
% varargout{1}: ld
% varargout{2}: ud
% 
%% PARSE THE INPUT DATA
narginchk(3,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('ld',[],@(x) (isnumeric(x)&numel(x)==numel(X))|isempty(x));
params.addParameter('ud',[],@(x) (isnumeric(x)&numel(x)==numel(X))|isempty(x));
params.addParameter('binscale','linear',@(x) x==1|x==0);
params.parse(varargin{:});

% Extract out variables from parsed input
LD=params.Results.ld;%lower error bar
UD=params.Results.ud;%upper error bar

% Check to make sure that X and Y are column arrays
if size(X,1)~=size(Y,1)
    error(['Binning aborted! Make sure that X and Y ',...
        'are column arrays with the same length.'])
end

% Check to make sure that if LD is not empty, UD has the same size
if ~isempty(LD)
   if numel(LD)~=numel(UD)&&numel(LD)==numel(X)
%        warning(['LD and UD have different sizes! ',...
%            'Make sure that ld is the same size as ud']);
   end
end

% Bin the data and take mean for each bin
if strcmp(params.Results.binscale,'linear')
    x=min(X):w:max(X);
elseif strcmp(params.Results.binscale,'log')
    x=10.^(log10(min(X)):w:log10(max(X)));
    x=x(find(~isnan(x)));
end
x=x';
x2=interp1(x,1:numel(x),X(:),'nearest');
ii1=find(~isnan(x2));
y=accumarray(x2(ii1),Y(ii1)',[numel(x) 1],@(x) mean(x));

if ~isempty(LD)
    ld=accumarray(x2(ii1),LD(ii1)',[numel(x) 1],@(x) mean(x));
    varargout{1}=ld;
end

if ~isempty(UD)
    ud=accumarray(x2(ii1),UD(ii1)',[numel(x) 1],@(x) mean(x));
    varargout{2}=ud;
end


