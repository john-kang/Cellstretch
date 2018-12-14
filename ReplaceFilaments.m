function [filNodeLUT_new]=ReplaceFilaments(intNodePositioning,fBreak,filNodeLUT,nodeXY,nIntNodes)
%ReplaceFilaments: Replace broken filaments (fBreak) with newly generated
%ones given node coordinates (nodeXY).


%11/07/28: repeats filaments!  Fix such that no filaments are repeated
%after replacement

if nargin==0
    clear all; close all;
    disp('Start ReplaceFilaments');
    intNodePositioning='r';
    fBreak=[1 3];
    filNodeLUT=[1 2; 2 3; 4 5]';
    nodeXY=[0 70.7106781186548,-70.7106781186547,-70.7106781186548,70.7106781186547;0 70.7106781186547,70.7106781186548,-70.7106781186547,-70.7106781186548];
    nIntNodes=1;
    PrintNetwork_XY_v2(nodeXY,filNodeLUT,nIntNodes);
end

nRemain=length(fBreak); %# filaments remaining to be generated
nNodes=size(nodeXY,2);

switch intNodePositioning
    case 'r' %internal nodes generated using randomly generated intersections
        filNodeLUT_break=filNodeLUT;
        filNodeLUT_break(:,fBreak)=[];
        filNodeLUT_new=GenerateFilaments(nRemain,nNodes,filNodeLUT_break);
        
    case 'i' %internal nodes generated using intersections
        %-
        %         nRemain=length(fBreak); %# filaments to be generated
        %         nPeriNodes=size(nodeXY,2)-nIntNodes; while nRemain>0
        %         %while there are still filaments to be generated...
        %             %filNodeLUT_periOnly =
        %             GenerateFilaments(nRemain,nPeriNodes);
        %
        %             [intXY,intxFils]=FindIntersections_v2(nIntNodes,periXY,filNodeLUT);
        %         end
        
        %just do the same thing as the random case for now
        filNodeLUT_break=filNodeLUT;
        filNodeLUT_break(:,fBreak)=[];
        filNodeLUT_new=GenerateFilaments(nRemain,nNodes,filNodeLUT_break);
      
end

if nargin==0
    figure;
    PrintNetwork_XY_v2(nodeXY,filNodeLUT_new,nIntNodes);
end

end