%This file contains functions:
%-SetupPeripheralNodes
%-MirrorQuadrants

%This file is called by:
%-GenerateNetwork

function [periPol,nPeriMobile,msg] = SetupPeripheralNodes(nPeriNodes,periNodePositioning,vertOffset,stretchMethod,cellRadius,doRotate)
%SetupPeripheralNodes: Setup the peripheral nodes using number of
%peripheral nodes (nPeriNodes), how they are positioned
%(periNodePositioning), what the vertical angle offset to determine mobile
%nodes is (vertOffset), the method of stretching (stretchMethod), the cell
%radius (cellRadius), whether to give the arrangement symmetry across
%both the X- and Y-axes (doMirror), and whether the peripheral nodes are
%rotated by 15 degrees (doRotate)

%Returns polar coordinates of peripheral nodes as well as how many of the
%vertical nodes are mobile (nPeriMobile)

%NOTE: the mobile node coordinates are listed first

%12/01/12: added doRotate parameter which rotates the peripheral nodes.
%Removed doMirror as an input parameter: is defaulted to True.

%12/01/06: added stretchMethod parameter which will determine, along with
%vertOffset, which mobile nodes are chosen. Added msg return parameter to
%return any error messages

%11/07/13: mirroring can be toggled now

switch nargin
    case 0
        disp('Start SetupPeripheralNodes');
        close all;
        nPeriNodes=8;
        periNodePositioning='u';
        cellRadius=100;
        vertOffset=45;
        stretchMethod='s';
        doRotate=false;
        doMirror=true;
    otherwise
        doMirror=true;
end

if doMirror && mod(nPeriNodes,4) ~= 0 && (periNodePositioning=='u'||periNodePositioning=='r') %uniform and random distributions may want bilateral symmetry
    msg='ERROR (SetupPeripheralNodes): nPeriNodes is not divisible by 4!';
    periPol=[]; nPeriMobile=[];
    disp(msg);
    return
end

%Place peripheral nodes at designated locations around circle
nPeriPerQuad=nPeriNodes/4;
switch(periNodePositioning)
    case('r')   %random positioning around circle
        if doMirror
            randFirstQuadTheta=rand(1,nPeriPerQuad)*90;   %assigns focal adhesions on cirle perimenter at random angles in quadrant #1
            [periDeg,msg]=MirrorQuadrants(randFirstQuadTheta);         %periPol(1,:) = RandomizePeripheralNodeAngles(nPeriNodes);
        else
            periDeg=rand(1,nPeriNodes)*360;
        end
        
    case('k')   %original positions used by Kathy
        periDeg=[80 90 100 260 270 280 0 30 60 120 150 180 210 240 300 330];%.*(pi/180);  %assigns specified angles for nPeriphNode # "focal adhesions"
        
    case('u')   %uniform position around circle
        
        if doMirror
            uniformFirstQuadTheta=(1/nPeriPerQuad *0.5 : 1/nPeriPerQuad : 1-1/nPeriPerQuad *0.5)*90;
            [periDeg,msg]=MirrorQuadrants(uniformFirstQuadTheta);
        else
            periDeg=(1/nPeriNodes/2:1/nPeriNodes:1-1/nPeriNodes/2)*360; %what is the difference???
        end
        
        if doRotate && mod(nPeriNodes,4)==0   %this is another way of having the uniform position with the nodes starting at angle 0 instead
            periDeg=0: 360/nPeriNodes :360;
            periDeg(end)=[];    %remove the last angle (360) which is the same as the first angle (0)
        end
end

switch(stretchMethod)
    case 's'
        %Fixed peripheral nodes are the ones that are stretched (apical side) or fixed
        %(basal side)
        fixedPeri=union(FindDataSubset(periDeg,90-vertOffset,90+vertOffset,1,'column data'),FindDataSubset(periDeg,270-vertOffset,270+vertOffset,1,'column data'));         %Find vertical nodes
        
        %Mobile peripheral nodes are the ones that are allowed to move to
        %relax stresses
        mobilePeri=1:nPeriNodes; mobilePeri(fixedPeri)=[];
    otherwise %currently 'p', 'i' methods
        mobilePeri=union(FindDataSubset(periDeg,90-vertOffset,90+vertOffset,1,'column data'),FindDataSubset(periDeg,270-vertOffset,270+vertOffset,1,'column data'));         %Find vertical nodes
        fixedPeri=1:nPeriNodes; fixedPeri(mobilePeri)=[];
end

nPeriMobile=length(mobilePeri);

%places mobile nodes at beginning and convert deg->radians 
periPol(1,:)=[periDeg(mobilePeri) periDeg(fixedPeri)]*pi/180;

%assigns radial location of focal adhesions on the perimeter
periPol(2,:) = cellRadius; 

if nargin==0    %for test case, graph
    [periXY(1,:) periXY(2,:)]=pol2cart(periPol(1,:), periPol(2,:));
    PrintNetwork_XY(periXY,[],0,cellRadius,[],[],{[periNodePositioning ' positioned peripheral nodes']; [stretchMethod ' stretch method']; ['mobile nodes (vertOffset ' num2str(vertOffset) ')  are cyan circles']},true,nPeriMobile);
%     hold on;
%     plot(periXY(1,1:nPeriMobile),periXY(2,1:nPeriMobile),'*','MarkerEdgeColor','c','MarkerFaceColor','c');  %plot mobile nodes as cyan circles
end

end

function [periTheta,msg]=MirrorQuadrants(firstQuadPeriNodeTheta)
%MirrorQuadrants: mirrors peripheral nodes in the 1st quadrant to all four
%quadrants
msg=[];
nPeriPerQuad = size(firstQuadPeriNodeTheta,2);
if size(firstQuadPeriNodeTheta,2)~=nPeriPerQuad
    msg='ERROR (SetupPeripheralNodes: MirrorQuadrants): size discrepancy between nPeriPerQuad and num peri nodes in first quadrant!';
    disp(msg)
    periTheta=[];
    return
end
periTheta=NaN(1,nPeriPerQuad*4);
periTheta(1:nPeriPerQuad) = firstQuadPeriNodeTheta;   %assigns focal adhesions on cirle perimenter at random angles in quadrant #1
periTheta(nPeriPerQuad+1:2*nPeriPerQuad) = 180-periTheta(1,1:nPeriPerQuad);  %quadrant #2
periTheta(2*nPeriPerQuad+1:3*nPeriPerQuad) = 180+periTheta(1,1:nPeriPerQuad); %quadrant #3
periTheta(3*nPeriPerQuad+1:4*nPeriPerQuad) = 360-periTheta(1,1:nPeriPerQuad); %quadrant #3
end

