function rsq=rsq_calc(x,y,p)
% Joshua Yeh
% Date: 19/05/27
% This fcn calculates the rsqaure value of a polynomial fit outputted by
% the polyfit MATLAB fcn.
%
%% INPUT VARIABLES
% x: xdata
% 
% y: ydata
% 
% p: the polynomial coefficients
% 
%%

% Calculate the mean in y
y_mean=mean(y);

% Calculate y-fit
y_fit=polyval(p,x);

% calc the regression sum of squares
SSR=sum((y_fit-y_mean).^2);

% calc the total sum of squares
SSTO=sum((y-y_mean).^2);

% calc r-squared values
rsq=SSR./SSTO;