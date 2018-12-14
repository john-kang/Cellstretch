function [Strain,strcFilAng] = CalcFilStrainsAndAngles(nodeXY,filNodeLUT,oldFilLen)
%CalculateFilStrains: calculates the filament strains after stretching
%   Calculate strain uses old filament lengths and current filament lengths.  Also calculates angle between the filaments.
filX(1,:) = nodeXY(1,filNodeLUT(1,:));      %\
filX(2,:) = nodeXY(1,filNodeLUT(2,:));      % \
filY(1,:) = nodeXY(2,filNodeLUT(1,:));      % -<2,nFilaments> with X,Y cds for stretched fil ends
filY(2,:) = nodeXY(2,filNodeLUT(2,:));      %/
dX = double(filX(2,:)-filX(1,:));           %-<1,nFilaments> calc delta <X,Y>
dY = double(filY(2,:)-filY(1,:));           %/
strcFilLen = double(sqrt(dX.^2 + dY.^2));    %<1,nFilaments> find new filaments length
strcFilAng = double(atan2(dY,dX));           %<1,nFilaments> find new filaments angle
delF = double(strcFilLen-oldFilLen);         %<1,nFilaments> oldFilLen was initialized first before stretching
Strain = double(delF./oldFilLen);            %<1,nFilaments> can be + or -
end

