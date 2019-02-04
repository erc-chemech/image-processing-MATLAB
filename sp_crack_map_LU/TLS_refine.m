function pp_new=TLS_refine(fit_obj,x,y)
%% DESCRIPTION
% This function refines a fit using a total least squares based on the pp
% structure (defined by a linear squares fit approach).
%
%% INPUT VARIABLES
% fit_obj: fit object from linear least squares fit
% x
pp=fit_obj.p;
fit_obj_struct=struct(fit_obj);
mypoly=@(x,p) polyval(p,x);%p order polynomial anon. fcn
c=pp.coefs;% get coefficients for each piece
pn=pp.pieces;%number of pieces
intervals=pp.breaks;%interval in which the polynomial was fitted
n=pp.order;%spline order
c_new2=[];
for dum=1:numel(intervals)-1
    c0=c(dum,:);%extract the polynomial coeffs for the local interval
    
    % extract the local points to fit the polynomial
    ii=x>=intervals(dum)&x<=intervals(dum+1);
    
    if sum(ii)<n
%         c_new2=cat(1,c_new2,c0);
        continue
    end
    x_int=x(ii);
    x_int_n=(x_int-fit_obj_struct.meanx)./fit_obj_struct.stdx;
    y_int=y(ii);
    
    [~,c_new,~]=numerFminS(mypoly,n,x_int_n,y_int,...
        'initial',c0,'res_f',10);
    pp.coefs(dum,:)=c_new;
%     c_new2=cat(1,c_new2,c_new);
end
pp_new=mkpp(intervals,c);
