function [fval, a,output]=numerFminS(fun, p, xdata, ydata,varargin)
%
% [Err, min_param]=numerFminS(fun, p, LBa, UBa, xdata, ydata)
%
% Orthogonal non-linear regression method in 2D for model 
% defined in file 'Model.m', which is an input as a 'fun'.
%     
% This function is based on optimization algorithm created
% by John D'Errico in 2006 and named as: 'fminsearchbnd'.
% 
% Input parameters:
%  - fun: name of file where is the model defined
%  - p: number of optimized parameters (same as in the file 'fun')
%  - LBa: lower bound vector or array for parameters 
%  - UBa: upper bound vector or array for parameters 
%  - xdata: input data block -- x: axis
%  - ydata: input data block -- y: axis
% 
% Return parameters:
%  - fval: error - sum of squared orthogonal distances 
%  - a: vector of model parameters  
%  - output: struct var returned from fminsearchbnd
%
% Authors:
% Ivo Petras (ivo.petras@tuke.sk)
% Tomas Skovranek (tomas.skovranek@tuke.sk)
% Dagmar Bednarova (dagmar.bednarova@tuke.sk)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was edited by Joshua Yeh (josh.yeh@espci.fr)
% 
% What changed:
% 
% - A weight factor was implemented in the fitting procedure. - Finding the
% minimum orthogonal distance was optimized by implementing the knnsearch
% MATLAB function and removing the for loop. A for loop is not needed to
% find the minimum distance between a data point and the cirve fit. This
% prevents repeated calls of knnsearch, which can slow down the fit
% significantly.
% 
% - Removed extraneous variables making the code more readable.
% 
% - The function outputs the fitting results and display the exitflag
% and the sum of the residuals.
% 
% - Placed the lower bound and upper bound inputs in to the varargin
% structure
% 
% - Added a "resolution factor", res_f, that controls the number of
% datapoints in the fitting curve. The higher this factor is, the more
% accurate the algorithm is in finding the point along the fitting curve
% that minimizes the orthogonal distance between the datapoint and the
% fitting curve. However, the higher the number, the longer the algorithm
% will take in determining the parameter fits.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date Created: 25/01/2009
% Date Edited: 6/01/2018
%
%
% Example:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [yM] = model(xm, a)
% yM = zeros(3,1);
% yM=a(1)*xm.^2+a(2)*xm+a(3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % measured data are stored in the file 'data.txt' 
% load data.txt;
% xdata=data(1,:);
% ydata=data(2,:);
% [ErrTLS,P]=numerFminS(@model,3,[-0.03 0.1 10.0], [-0.15 0.9 11], xdata, ydata)
% YhatTLS=polyval(P2(1:3),xdata);
% plot(xdata, ydata, '*');
% hold on
% plot(xdata,YhatTLS,'k');
%

% Parse the varargin
narginchk(4,inf); %Check the number of input values
params=inputParser;
params.CaseSensitive=false;
params.addParameter('wx',[]);
params.addParameter('wy',[]);
params.addParameter('initial',zeros(1,p),@(x) length(x)==p);
params.addParameter('LB',ones(1,p).*-inf,@(x) length(x)==p);
params.addParameter('UB',ones(1,p).*inf,@(x) length(x)==p);
params.addParameter('res_f',4,@(x) isnumeric(x));
params.parse(varargin{:});

wx=params.Results.wx;%weight factors in x
wy=params.Results.wy;%weight factors in y
a0=params.Results.initial;%initial parameter guess
LBa=params.Results.LB;%lower bound
UBa=params.Results.UB;%upper bound
res_f=params.Results.res_f;%resolution factor of the fitting curve

if isempty(wx)%check to see if weight factor was set
    wx=ones(numel(xdata),1);
end

if isempty(wy)%check to see if weight factor was set
    wy=ones(numel(ydata),1);
end

if ~(exist('fminsearchbnd', 'file') == 2)
    P = requireFEXpackage(8277);  
    % fminsearchbnd is part of 8277 at MathWorks.com
end 

warning off all
options = optimset('MaxIter',1e+4,'MaxFunEvals',1e+4,'TolX',1e-6,'TolFun',1e-6); 

Mpoints=[xdata,ydata];
xx=linspace(min(xdata), max(xdata), numel(xdata)*res_f)';
[a,fval,exitflag,output] = fminsearchbnd(@calculation,a0,LBa,UBa,options);

% Display the results
disp(['exitflag: ',num2str(exitflag),'   residual sum: ',num2str(fval)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function for evaluating the model with optimized parameters
    function [sum1] = calculation (a)

        [yy]=fun(xx, a);
        points=[xx yy];

        [ii,~]=knnsearch(points,Mpoints);
        w_min_dist=sqrt(wx.*(xx(ii)-Mpoints(:,1)).^2+...
            wy.*(yy(ii)-Mpoints(:,2)).^2);%weighted minimum distance
        sum1=sum(w_min_dist);%sum of residuals in both x and y
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%