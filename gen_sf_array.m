function out=gen_sf_array(sf,varargin)
% Author: Joshua Yeh
% Date created: 18/03/19
% DESCRIPTION: This function takes a surface fit object and generates a
% correction square image array based on a user-defined resolution size of
% the image array.
%% INPUT
% SF: surface fit object
% varargin
%      'type': The square image resolution (256 or 512 or 1024 or 2048)
%      default is 512, NOTE that this assumes that the surface fit object
%      is calculated based on a 2048 by 2048 object
%% OUPUT
% out: the output surface map
%%
% Parse the input variables
narginchk(1,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('L',1:4:2048);
params.parse(varargin{:});

% Extract out values from parsed input
L=params.Results.L;

% Create the surface map
disp('This can take awhile, ~20 s');
out=feval(sf,repmat(L,numel(L),1),repmat(L,1,numel(L)));
disp('Finished generating surface array');
