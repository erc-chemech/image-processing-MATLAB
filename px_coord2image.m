function out=px_coord2image(rows,cols,zdata,I_size)
% Author: Joshua Yeh
% Date created: 18/06/04
%% DESCRIPTION
% This function accepts the corresponding pixel row and columns and maps
% it into an image array defined by I_size;
%
%% INPUT VARIABLES
% row: row pixel coordinate
% col: column pixel coordinate
% I_size: size of the image array
% 
%% OUTPUT VARIABLES
% out: the image array with elements containing values in zdata with row
% and column coordinates defined by its corresponding pixel coordinate
% 
%%

% Make everything to col array
[rows,cols,zdata]=prepareSurfaceData(rows,cols,zdata);

out=nan(max(rows),max(cols));% preallocate image array

for dum=1:numel(rows)
    r=rows(dum);%px row
    c=cols(dum);%px column
    z=zdata(dum);%zdata at r, c
    out(r,c)=z;%store the zdata @ corresponding r, c
end