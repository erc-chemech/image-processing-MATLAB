function center_axes(ax,varargin)
% Author: Joshua Yeh
% Date created: 2018/05/11
% 
%% DESCRIPTION
% This function accepts an axes handle and recenters the axes in the figure
% window.
% 
%% INPUT VARIABLES
% ax: axes handles
%
% varargin: <field name>, '<key>'
% 
% 'margins': additional margins in units of pixels
% 
%%

    %  Parse input variables
    if nargin<1
        ax=gca;
    end
    
    params=inputParser;
    params.CaseSensitive=false;
    params.addParameter('margins',25,@(x) isnumeric(x));
    params.addParameter('c_dy',5,@(x) isnumeric(x));
    params.parse(varargin{:});
    
    % Extract out values from parse input
    c=params.Results.margins;
    c_dy=params.Results.c_dy;

    if isa(ax,'matlab.graphics.axis.Axes')
        ca(ax,c,c_dy);
    elseif isa(ax,'matlab.ui.Figure')
        axs=findall(ax,'type','axes');
        for dum=1:numel(axs)
            ca(axs(dum),c);
        end
    end

    
end

function ca(ax,c,c_dy)
    fig=ax.Parent;%figure handle
    fig.Units='pixels';
    ax.Units='pixels';
    fip=fig.InnerPosition;%figure inner position
    ti = ax.TightInset; %margins of text labels
    outerpos = ax.OuterPosition;%outer position of the axes
    c_dx=0;%additional x offset needed if colorbar exists
%     c_dy=5;%additional y offset needed if colorbar exists
    c_dw=0;
    c_dh=0;
    
    % Check to see if there is a colorbar associated with the axes
    if ~isempty(ax.Colorbar)
        cb=ax.Colorbar;%colorbar handles
        cL=cb.YLabel;%colorbar label
        cb.Units='pixels';
        cL.Units='pixels';
        mh=max([ti(2) ti(4)]);
        mv=max([ti(1) ti(3)]);
        
        % Determine additional offset needed so that the
        % colorbar is included
        switch cb.Location
            case 'southoutside'
                    
                if cL.Extent(2)<-mh
                    
                    c_dy=ax.Position(2)-ti(2)-(cb.Position(2)+cL.Extent(2))+5;
                    c_dh=c_dy;
                    
                elseif cL.Extent(2)>-mh&&strcmp(cb.YAxisLocation,'bottom')
                    
                    c_dy=ax.Position(2)-ti(2)-cb.Position(2)+mh;
                    c_dh=c_dy;
                   
                else
                    c_dy=ax.Position(2)-ti(2)-cb.Position(2)+5;
                    c_dh=c_dy;
                    
                end
                
            case 'eastoutside'
                
                if cb.Position(1)+cL.Extent(1)+cL.Extent(3)>...
                        cb.Position(1)+mv
                    
                    c_dw=cb.Position(1)+cL.Extent(1)+cL.Extent(3)-...
                        (ax.Position(1)+ax.Position(3)+ti(3))+5;
                    
                elseif cb.Position(1)+cL.Extent(1)+cL.Extent(3)<...
                        cb.Position(1)+mv&&...
                        strcmp(cb.YAxisLocation,'right')
                    
                    c_dw=cb.Position(1)+mv-...
                        (ax.Position(1)+ax.Position(3)+ti(3));
                    
                else
                    
                    c_dw=cb.Position(1)+mv-...
                        (ax.Position(1)+ax.Position(3)+ti(3));
                    
                end
                
            case 'northoutside'
                
                if cb.Position(2)+cL.Extent(2)+cL.Extent(4)>...
                        cb.Position(2)+mh
                    
                    c_dh=cb.Position(2)+cL.Extent(2)+cL.Extent(4)-...
                        (ax.Position(2)+ax.Position(4)+ti(4))+5;
                    
                elseif cb.Position(2)+cL.Extent(2)+cL.Extent(4)<...
                        cb.Position(2)+mh&&strcmp(cb.YAxisLocation,'top')
                    
                    c_dh=cb.Position(2)+mh-...
                        (ax.Position(2)+ax.Position(4)+ti(4))+5;
                    
                else
                    
                    c_dh=cb.Position(2)+cb.Position(4)-...
                        (ax.Position(2)+ax.Position(4)+ti(4))+5;
                    
                end
                
            case 'westoutside'
                
                if cb.Position(1)+cL.Extent(1)>cb.Position(1)-mh
                    
                    c_dx=ax.Position(1)-ti(1)-(cb.Position(1)+cL.Extent(1))+5;
                    c_dw=c_dx+mh/2;
                    
                elseif cb.Position(1)+cL.Extent(1)<cb.Position(1)-mh&&...
                        strcmp(cb.YAxisLocation,'left')
                    
                    c_dx=ax.Position(1)-ti(1)-(cb.Position(1)-mh)+5;
                    c_dw=c_dx+mh/2;
                    
                else
                    
                    c_dx=ax.Position(1)-ti(1)-cb.Position(1)+5;
                    c_dw=c_dx+mh/2;
                    
                end
                
        end
    end
    
    
    
    
    % find target x position in pxs
    targetx=ti(1)+c_dx+c;
    
    % find target y position in pxs
    targety=ti(2)+c_dy+c;
    
    % find target width
    targetw=fip(3)-ti(1)-ti(3)-c_dw-2*c;
    
    % find target height
    targeth=fip(4)-ti(2)-ti(4)-c_dh-2*c;
    
    ax.Position = [targetx targety targetw targeth];
    
    ax.Units='normalized';
    try cb.Units='normalized';end
end