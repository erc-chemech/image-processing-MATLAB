function ax=alphaColorbar_stacked(ref_ax,max1,cbs)
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

levels=0:max1;
levels=cat(1,levels,levels,levels);
levels_n=levels/max1;
ap=discretize(levels_n,ax.Parent.Alphamap);

[X,Y]=meshgrid(levels(1,:),[0 0.5 1]);
for dum=1:numel(cbs)
    surf(ax,X,Y+dum-1,ones(size(X)),'alphadatamapping','direct',...
    'facealpha','flat','alphadata',ap,...
    'facecolor',cbs{dum},'cdata',levels);
    if dum<numel(cbs)
        plot(ax,[0 max1],[dum dum],'k-','linewidth',1);
    end
end



set(ax,'box','on','xlim',[0 max1],'ylim',[0 numel(cbs)],...
    'fontname',ref_ax.FontName,'fontsize',ref_ax.FontSize,...
    'xgrid','off','ygrid','off','zgrid','off','ytick',[],...
    'view',[0 90],'ycolor','k');
xlabel(ax,'Counts','fontsize',ax.FontSize);
