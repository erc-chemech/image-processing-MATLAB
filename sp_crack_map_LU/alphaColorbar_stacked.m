function ax=alphaColorbar_stacked(ref_ax,max1)
%This funcion creates a custom alpha colorbar for Figure 7
% 
%% INPUT VARIABLES
% ref_ax: reference axes
% 
% max1: maximum counts
%
%% OUTPUT VARIABLES
%  ax: axes of the custom stacked colorbar
% 

ax=axes;
ax.Parent=ref_ax.Parent;
hold(ax,'on');
ax.Position=[ref_ax.Position(1) 0.15 ref_ax.Position(3) 0.15];
ax.Parent.Alphamap=linspace(0,1,256);

max2=max1*1e1;

% loading colorbar
levels_l=1:max2;
levels_l=cat(1,levels_l,levels_l,levels_l);
levels_n_l=levels_l/max2;
ap_l=discretize(levels_n_l,ax.Parent.Alphamap);
[X,Y]=meshgrid(levels_l(1,:),[0 0.5 1]);
surf(ax,X,Y+2,ones(size(X)),'alphadatamapping','direct',...
    'facealpha','flat','alphadata',ap_l,...
    'facecolor','b','cdata',levels_l);

% unloading colorbar
levels_u=1:max1;
levels_u=cat(1,levels_u,levels_u,levels_u);
levels_n_u=levels_u/max1;
ap_u=discretize(levels_n_u,ax.Parent.Alphamap);
[X,Y]=meshgrid(levels_u(1,:),[0 0.5 1]);
surf(ax,X,Y+1,ones(size(X)),'alphadatamapping','direct',...
    'facealpha','flat','alphadata',ap_u,...
    'facecolor','m','cdata',levels_u);

[Xx,Yy]=meshgrid(max1:max2,[0 0.5 1]);
surf(ax,Xx,Yy+1,ones(size(Xx)),'facecolor','m');

% "other" colorbar
surf(ax,X,Y,ones(size(X)),'alphadatamapping','direct',...
    'facealpha','flat','alphadata',ap_u,...
    'facecolor','k','cdata',levels_u);

[Xx,Yy]=meshgrid(max1:max2,[0 0.5 1]);
surf(ax,Xx,Yy,ones(size(Xx)),'facecolor','k');

set(ax,'box','off','xlim',[0 max2],'ylim',[0 3],...
    'fontname',ref_ax.FontName,'fontsize',ref_ax.FontSize,...
    'xgrid','off','ygrid','off','zgrid','off','ytick',[],...
    'view',[0 90],'ycolor','none','xscale','log','tickdir','out');
xlabel(ax,'Counts','fontsize',ax.FontSize);
