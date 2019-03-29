function out=my_honeycomb(ax,xx,yy,varargin)

%% Parse inpurt variables
narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('res','auto',@(x) isnumeric(x));
params.addParameter('max',100,@(x) isnumeric(x));
params.addParameter('color','k');
params.addParameter('EdgeColor','none');
params.addParameter('dn',[]);
params.addParameter('xEdges',[],@(x) isnumeric(x));
params.addParameter('yEdges',[],@(x) isnumeric(x));
params.parse(varargin{:});

% Extract out balues from parsed input
res=params.Results.res;
max1=params.Results.max;
max1=round(max1);
color=params.Results.color;
edgecolor=params.Results.EdgeColor;
dn=params.Results.dn;

ax.Parent.Alphamap=linspace(0,1,256);

if ~ischar(res)
    if isempty(res)
        % determine step resolution
        xe=linspace(ax.XLim(1),ax.XLim(2),100);
        ye=linspace(ax.YLim(1),ax.YLim(2),100);
        res=mean([xe(2)-xe(1) ye(2)-ye(1)]);
    end
    if isempty(params.Results.xEdges)
        xEdges=ax.XLim(1):res:ax.XLim(2);%x edges
    else
        xEdges=params.Results.xEdges(1):res:params.Results.xEdges(2);
    end
    if isempty(params.Results.yEdges)
        yEdges=ax.YLim(1):res:ax.YLim(2);%y edges
    else
        yEdges=params.Results.yEdges(1):res:params.Results.yEdges(2);
    end
    out=honeycomb(xx,yy,'ax',ax,'xEdges',xEdges','yEdges',yEdges');
elseif strcmp(res,'auto')
    out=honeycomb(xx,yy,'ax',ax);
end

N=out.CData;
N_norm=N./max1;%normalize 2d count array
N_norm(N_norm>1)=1;
ap=discretize(N_norm,ax.Parent.Alphamap);
out.DisplayName=dn;

set(out,'edgecolor',edgecolor,'facecolor',color,'AlphaDataMapping','direct',...
    'facealpha','flat','facevertexalphadata',ap,'userdata',max1);

% keyboard