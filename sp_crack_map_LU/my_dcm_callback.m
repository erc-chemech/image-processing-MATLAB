function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

ax=event_obj.Target.Parent;
line_handles=findall(ax,'type','line');

xdatas=[];
for dum=1:numel(line_handles)
    xdata=line_handles(dum).XData;
    xdatas=cat(1,xdatas,xdata');
end
xdatas=unique(xdatas);
xdatas=sort(xdatas);
ii=find(xdatas==pos(1));
output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)],...
    ['index: ',num2str(ii)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
