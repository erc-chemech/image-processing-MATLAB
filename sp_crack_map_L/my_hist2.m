function surf1=my_hist2(ax,xx,yy,varargin)
%% Parse inpurt variables
narginchk(2,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('res',100,@(x) isnumeric(x));
params.addParameter('max',100,@(x) isnumeric(x));
params.addParameter('color','k');
params.addParameter('plot','on',@(x) ischar(x));
params.addParameter('EdgeColor','none');
params.parse(varargin{:});

% Extract out balues from parsed input
res=params.Results.res;
max1=params.Results.max;
max1=round(max1);
color=params.Results.color;
flag=params.Results.plot;
edgecolor=params.Results.EdgeColor;

ax.Parent.Alphamap=linspace(0,1,256);

x=linspace(ax.XLim(1),ax.XLim(2),res);%x edges
y=linspace(ax.YLim(1),ax.YLim(2),res);%y edges
[N,~,~]=histcounts2(xx(:),yy(:),x,y);
N=N';

N_norm=N./max1;%normalize 2d count array
N_norm(N_norm>1)=1;
ap=discretize(N_norm,ax.Parent.Alphamap);
N(N==0)=nan;% remove 0 components

x=x(1:end-1);%associate columns to correct BCC values
y=y(1:end-1)+abs(y(2)-y(1));%associate rows to correct RCC values
[Y,X]=meshgrid(x,y);

if strcmp(flag,'on')
    surf1=surf(ax,Y,X,N,'alphadatamapping','direct',...
        'facecolor',color,'facealpha','flat',...
        'alphadata',ap,'userdata',max1,'edgecolor',edgecolor);
else
    surf1.N=N;
    surf1.X=Y;
    surf1.Y=X;
end