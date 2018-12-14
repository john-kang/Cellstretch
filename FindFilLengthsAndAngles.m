function [ filLens, filAng ] = FindFilLengthsAndAngles(nodeXY,filNodeLUT)
%FindFilLengthsAndAngles:
%   Find the lengths and angles (using atan2) of filaments using node XY coords and the lookup
%   table for nodes and filaments

%11/09/16:

if nargin==0
    nodeXY=[0 0; sqrt(pi)/2 sqrt(pi)/2]';
    filNodeLUT=[1 2]';
end

filX(1,:)=nodeXY(1,filNodeLUT(1,:));    %x-coord of origin nodes
filX(2,:)=nodeXY(1,filNodeLUT(2,:));    %x-coord of destination nodes
filY(1,:)=nodeXY(2,filNodeLUT(1,:));    %y-coord of origin nodes
filY(2,:)=nodeXY(2,filNodeLUT(2,:));    %y-coord of destination node

dX = double(filX(2,:)-filX(1,:));       %change in x-coordinates for filament
dY = double(filY(2,:)-filY(1,:));       %change in y-coordinates for filament 

filLens = double(sqrt(dX.^2 + dY.^2));  %length of filament

filAng = double(atan2(dY,dX));           %<1,nFilaments> find new filaments angle

end