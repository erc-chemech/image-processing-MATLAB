
function [load,unload]=parse_scl(time,stress,varargin)
% Joshua Yeh
% Date: 190521
%% DESCRIPTION
% This fcn parses step cyclic loading profile into loading & unloading.
% 
%% VARIABLE INPUTS
% varargin (variable argument input)
% 'fieldname',<fieldvalue>
%
% 'show': show the results of the parsing analysis
% 
% 'include_load_end': include the last datapoint as part of peak stresses
%
% 'include_unload_end': unclude the last datapt as part of min. or 0 stress
% 
% 'minpeakpriminence': defines the min peak prominence (see findpeaks fcn)
% 
%% OUTPUT VARIABLES
% 'load': struct var containing indices corres. to loading
% 
% 'unload': struct va containing indices corres. to unloading
% 
%% Parse inputs
narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('show','off',@(x) strcmp(x,'on')|strcmp(x,'off'));
params.addParameter('include_load_end',1,@(x) isnumeric(x));
params.addParameter('include_unload_end',0,@(x) isnumeric(x));
params.addParameter('minpeakprominence',0.4,@(x) isnumeric(x));
params.parse(varargin{:});

%% find peaks
[~,ii]=findpeaks(stress,...
    'minpeakprominence',params.Results.minpeakprominence);
ii=[1;ii];% include first point
if params.Results.include_load_end==1
    ii=[ii;numel(stress)];%include last pt
end

%% find valleys
[~,ii2]=findpeaks(-stress+100,...
    'minpeakprominence',params.Results.minpeakprominence);
ii2=[1;ii2];% include first point
if params.Results.include_unload_end==1
    ii2=[ii2;numel(stress)];%include last pt
end

%% Parse loading
for dum=1:length(ii)-1
    stress_r=stress(ii2(dum):ii(dum+1));% current stress range
    load{dum,1}=[find(stress_r>stress(ii(dum)))+ii2(dum)-1];
end
disp('Identified loading portions!');

%% Parse unloading
for dum=1:length(ii2)-1
    stress_r=ii(dum+1):ii2(dum+1);
    unload{dum,1}=stress_r';
end
disp('Identified unloading portions!');

%% plot the parsing results
if strcmp(params.Results.show,'on')
    pf1=my_fig(99,{[1 1 1]});
    plot(pf1.s1,time,stress,'k.');%plot raw data

    for dum=1:length(load)%plot loading portions
        plot(pf1.s1,time(load{dum}),stress(load{dum}),'r.')
    end

    for dum=1:length(unload)%plot unloading portions
        plot(pf1.s1,time(unload{dum}),stress(unload{dum}),'b.');
    end

    xylabels(pf1.s1,'time (s)','nominal stress (MPa)');
    center_axes(pf1.s1);
end
