function [urx,ury,GCC_0_stress,residuals]=...
    gen_unload_array(unload_coef_fit1,unload_coef_fit2,unload_coef_fit3,...
    unload_coef_fit4,load_fit_TN135,GCC_0_stress)
%% DESCRIPTION
% This function generates the 2d unloading array that can be used to
% associated unloading chromatic poits to the maximum load history.
% 
%% INPUT VARIABLES
% unload_coef_fit1: fit object for the first 3rd order polynomial coeff.
% used to trace the unloading curve
% 
% unload_coef_fit2: fit object for the 2nd 3rd order polynomial coeff.
% used to trace the unloading curve
% 
% unload_coef_fit3: fit object for the 3rd 3rd order polynomial coeff.
% used to trace the unloading curve
% 
% unload_coef_fit4: fit object for the 4th 3rd order polynomial coeff.
% used to trace the unloading curve
% 
% load_fit_TN135
% 
%% OUTPUT VARIABLES
% urx: GCC 2D array
% 
% ury: TCC 2D array
% 
%%

% Conditions for performing algorithm
GCC_X0=1e-4;%initial guess for intersection
GCC_X=GCC_X0;
n1=5e4;% # pts in GCC of unloading curve
n2=5e4;% # pts in TCC of loading curve
que=250;%query points for calc abs dist between unloading and loading points
tol2=1e-6;% tolerance level req. to accept unloading curve

% Create loading chromatic curve
load_TCC=linspace(0,0.1297,n2)';
load_GCC=load_fit_TN135(load_TCC);
[load_TCC,ii]=sort(load_TCC);
load_GCC=load_GCC(ii);

% Preallocate 2d arrays
urx=nan(1e3,numel(GCC_0_stress));%GCC value of unloading array
ury=nan(1e3,numel(GCC_0_stress));%TCC value of unloading array
residuals=nan(numel(GCC_0_stress),1);%final "residual" of all unloading curves
residual_m=nan;%memory of last "residual" calculation

flag2=1;% toggle condition for outer while loop

for dum=1:numel(GCC_0_stress)
    
    if GCC_X>=5e-4
    else
        GCC_X=GCC_X0;
    end
    
    % Fit parameters that describe the polynomial fit of the unloading
    % curve
    p1=unload_coef_fit1(GCC_0_stress(dum));
    p2=unload_coef_fit2(GCC_0_stress(dum));
    p3=unload_coef_fit3(GCC_0_stress(dum));
    p4=unload_coef_fit4(GCC_0_stress(dum));
    params=[p1 p2 p3 p4];
    
    % Determine unloading curve based on initial guess
    temp_GCC=linspace(GCC_0_stress(dum),GCC_X,n1)';
    temp_GCC_n=(temp_GCC-mean(temp_GCC))./std(temp_GCC);% normalize xdata
    temp_TCC=polyval(params,temp_GCC_n);
    
    while flag2==1
        
        % find GCC values within a defined tolerance of load_GCC_q in the
        % unloading curve
        GCC_q=temp_GCC(end-que:2:end);
        TCC_q=temp_TCC(end-que:2:end);
        
        % Estimate the intersection point by calc. nearest neighbor
        [ii1,d]=knnsearch([GCC_q TCC_q],[load_GCC load_TCC]);
        [~,ii2]=min(d);
        load_GCC_q=load_GCC(ii2);
        
        % Determine "residual" based on the diff b/w max unloading GCC and
        % GCC at intersection
        GCC_q=sort(GCC_q);
        residual=abs(GCC_q(end)-load_GCC_q);
                
        if residual<tol2%cond for accepting unloading curve parameters
            
            disp(['fitted unloading curve #',num2str(dum),...
                'diff/tol2: ',num2str(residual/tol2)]);
            residuals(dum)=residual;%store final residual of unloading curve
            break
            
        else
            
            if residual==residual_m%cond for det no improvement in iteration
                
                fprintf([num2str(residual/tol2),...
                    '\t No improvement made with further iterations.\n']);
                residuals(dum)=residual;%store final residual of unloading curve
                break
                
            elseif residual>residual_m%cond for case where iteration is worse
                
                GCC_X=GCC_X_m;
                residuals(dum)=residual_m;
                break
                
            else%cond to continue optimizing unloading curve

                % redefine estimate of intersection point GCC_X
                GCC_X_m=GCC_X;% first, store memory of the "old" GCC_X
                GCC_X=GCC_q(ii1(ii2));

                % Update the unloading curve
                temp_GCC=linspace(GCC_0_stress(dum),GCC_X,n1)';
                temp_GCC_n=(temp_GCC-mean(temp_GCC))./std(temp_GCC);
                temp_TCC=polyval(params,temp_GCC_n);
                
            end
        
            residual_m=residual;%remember residual value
            
        end
        
    end

    % Update the unloading curve
    temp_GCC1=linspace(GCC_0_stress(dum),GCC_X,1e3)';%reduce the resolution
    temp_GCC_n=(temp_GCC1-nanmean(temp_GCC))./nanstd(temp_GCC);%normalize
    
    urx(:,dum)=temp_GCC1;
    ury(:,dum)=polyval(params,temp_GCC_n);
    
    % reset residual memory
    residual_m=nan;
    
end
