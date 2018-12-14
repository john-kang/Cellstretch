%This file is called by:
%-GenerateNetwork

function [periPol,nPeriMobile,msg] = SetupPeripheralNodes_Square(nPeriNodes,periNodePositioning,vertOffset,stretchMethod,cellRadius,doRotate)
%SetupPeripheralNodes_Square: 

switch nargin
    case 0
        disp('Start SetupPeripheralNodes');
        close all;
        nPeriNodes=60;
        periNodePositioning='u';
        cellRadius=100;
        vertOffset=45;
        stretchMethod='s';
        doRotate=false;
        doMirror=true;
    otherwise
        doMirror=true;
end
        
%% for flat network per SR reviewr
periNodePositioning = 's';  %square
nPeriSide=nPeriNodes/4;
spacing = linspace(-cellRadius, cellRadius, nPeriSide+2);
spacing(1) = [];    %remove ends
spacing(end) = [];   %remove ends

%arrange peripheral nodes in such a way: example of 8 perinodes in a flat
%network.  Top and bottom nodes * are fixed.  Left and right nodes o are mobile
%  ---*---*---
% |           |
% o           o
% |           |
% o           o
% |           |
%  ---*---*---

%mobile peripheral nodes first
periXY(1,1:nPeriSide) = -cellRadius;   %X-axis, LEFT side
periXY(2,1:nPeriSide) = spacing;   %Y-axis, LEFT side

periXY(1,nPeriSide+1:2*nPeriSide) = cellRadius;   %X-axis, RIGHT side
periXY(2,nPeriSide+1:2*nPeriSide) = spacing;   %Y-axis, RIGHT side

%fixed peripheral nodes second
periXY(1,2*nPeriSide+1:3*nPeriSide) = spacing;   %X-axis, TOP side
periXY(2,2*nPeriSide+1:3*nPeriSide) = cellRadius;   %Y-axis, TOP side

periXY(1,3*nPeriSide+1:4*nPeriSide) = spacing;   %X-axis, BOT side
periXY(2,3*nPeriSide+1:4*nPeriSide) = -cellRadius;   %Y-axis, BOT side

plot(periXY(1,:), periXY(2,:), '*');

[periPol(1,:), periPol(2,:)] = cart2pol(periXY(1,:),periXY(2,:));

nPeriMobile = nPeriSide*2;  %LEFT and RIGHT sides are mobile

%% for test case, graph
% if nargin==0    
    [periXY(1,:), periXY(2,:)]=pol2cart(periPol(1,:), periPol(2,:));
    PrintNetwork_XY(periXY,[],0,cellRadius,[],[],{[periNodePositioning ' positioned peripheral nodes']; [stretchMethod ' stretch method']; ['mobile nodes (vertOffset ' num2str(vertOffset) ')  are cyan circles']},true,nPeriMobile);
%     hold on;
%     plot(periXY(1,1:nPeriMobile),periXY(2,1:nPeriMobile),'*','MarkerEdgeColor','c','MarkerFaceColor','c');  %plot mobile nodes as cyan circles
% end

end
