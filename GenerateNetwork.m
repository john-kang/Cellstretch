%This file contains:
% -GenerateNetwork
% -Appendix (alternate code)

%This file called by cellstretch_john_execute

function [intPol,periPol,nPeriMobile,filNodeLUT]=GenerateNetwork(nIntNodes,nPeriNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,nFilaments,cellRadius,intNodeNoise,doRotate)
%GenerateNetwork: generates the actin network model given the basic
%parameters.  Model returned is defined by:
% 1,2) polar crds of the internal and peripheral nodes
% 3) # of mobile% peripheral nodes (given first before the fixed ones)
% 4) filaments, which are defined as connections between nodes


%14/12/22: added SetupPeripheralNodes_Square to take into acount a square
%configuration for the filament network

%12/01/12: updateed the filament generator to take into account the
%"rotated" uniform peripheral network that startes the peripheral nodes at
%angle 0.  Previously, the max number of filaments is not being generated
%(see the "Final Network" figure displayed).

%12/01/06: takes into account stretchMethod for determining whether the
%nodes described in vertOffset are to remain mobile, be stretched, or be
%stationary.  In methods besides the Stossel 2011 method, vertOffset
%describes the peripheral nodes that are to remain mobile.  With the
%Stossel method, vertOffset describes the nodes that are to be stretched
%(apical side) and the nodes that are to be held stationary (basal side).

%11/12/30: takes into account Stossel 2011 method for int node positioning.

%11/07/25: both random and intersect cases appear to generate correct
%networks.  For the intersect case, the algorithm will return error if the
%desired # intx (nIntNodes) cannot be obtained, which is dependent on how
%many filaments are initially generated.  Specifically, the nIntNodes
%cannot be greater than (total # filaments - 2*nIntNodes) choose 2. Other
%input error checks are also performed

%11/07/07: works but nFilaments is inconsistent with nIntNodes for the
%intersection case currently

%11/06/29: initial version

switch nargin
    case 0
        close all;
        clear all;
        disp('***Start GenerateNetwork***')
        
        intNodePositioning='s';
        periNodePositioning='u';    %"square" per SR reviewer
        vertOffset=45;
        stretchMethod='s';
        cellRadius=100;
        doRotate=true;
        
        nPeriNodes=20;
        
        if intNodePositioning=='s' && periNodePositioning=='u'
            %Calculate appropriate # internal nodes for a Stossel network
            nIntNodes=0;
            for ii=1:nPeriNodes/4 %for one hemisphere
                nIntNodes=nIntNodes+2*(ii-1);
            end
            nIntNodes=2*nIntNodes; %for other hemisphere
            
            if doRotate
                nIntNodes=nIntNodes+1;  %there is extra node possible when peripheral nodes are rotated
            end
            
            %Special exception for a 4 peripheral node network
            if nPeriNodes==4
                nIntNodes=1;
                nFilaments=8;
            end
        elseif intNodePositioning=='s' && periNodePositioning=='s'
            nIntNodes = (nPeriNodes/4)^2;
        end
        
        connectDensity = 0.5;
        
        nFilaments = CalcNumFilaments(intNodePositioning,nIntNodes,nPeriNodes,connectDensity,doRotate, periNodePositioning); %45;
        
        intNodeNoise = 0; %no noise
        
        dbstop if error
    otherwise
        %         doMirror=true;
end

er=1e-8; %number that approximates zero to avoid rounding errors

%SETUP PERIPHERAL NODES
%NOTE: the first "nPeriMobile" peripheral nodes are "mobile" like internal nodes


if periNodePositioning=='s';
    [periPol,nPeriMobile]=SetupPeripheralNodes_Square(nPeriNodes,periNodePositioning,vertOffset,stretchMethod,cellRadius,doRotate);
else
    [periPol,nPeriMobile]=SetupPeripheralNodes(nPeriNodes,periNodePositioning,vertOffset,stretchMethod,cellRadius,doRotate);
end

%SETUP INTERNAL NODES AND FILAMENTS
switch intNodePositioning
    
    case{'r'}   %randomly generate nodes on the perimeter
        %         [intPol]=SetupInternalNodes(nIntNodes,intNodePositioning,cellRadius);
        %         %initial location of internal nodes in polar coords,
        %         =nodes
        filNodeLUT = GenerateRandomFilaments(nFilaments,nIntNodes+nPeriNodes);  %generate filaments randomly using all nodes
        intPol(1,:) = rand(1,nIntNodes)*2*pi;   %random assignment of angles to nodes
        intPol(2,:) = rand(1,nIntNodes)*cellRadius;  %random assignment of radius
        
    case{'i','s'}   %use the intersection method or use the Stossel method (will diverge later)
        
        nPeriFilaments=nFilaments-2*nIntNodes;  %# filaments initally generated from peri node to peri node.  nFilaments will be the final # filaments. 2*nIntNodes is subtracted because more filaments will be generated during the intersection process
        
        %ERROR CHECKS
        if nPeriFilaments > nchoosek(nPeriNodes,2)
            disp(['ERROR (GenerateNetwork): desired # filaments (',int2str(nPeriFilaments),') not possible.  Max # is ',int2str(nchoosek(nPeriNodes,2))]);
            return
%         elseif nPeriFilaments < 2
%             disp(['ERROR (GenerateNetwork): initial # filaments (',int2str(nPeriFilaments),') generates a trivial network.  Increase #.']);
%             return
        elseif nIntNodes > nchoosek(nPeriFilaments,2)
            disp(['ERROR (GenerateNetwork): desired # internal nodes (',int2str(nIntNodes),') not possible.  Max # is ',int2str(nchoosek(nPeriFilaments,2))]);
            return
        end
        
        [periXY(1,:), periXY(2,:)]=pol2cart(periPol(1,:),periPol(2,:)); %convert peripheral node positions to polar crds
        
        switch intNodePositioning
            case 'i'
                nIntx=-1;
                trial=0;
                while nIntx<nIntNodes %while the # int nodes not sufficient, loop
                    %Note: passing Error Checks guarantees that desired # int nodes will
                    %be found
                    
                    filNodeLUT_periOnly = GenerateRandomFilaments(nPeriFilaments,nPeriNodes);  %generate filaments randomly using only peripheral nodes.
                    
                    [intXY,intxFils] = FindIntersections_Schwarz(nIntNodes,periXY,filNodeLUT_periOnly);    %generate nIntNodes number of internal nodes at intersections of filaments
                    
                    nIntx=size(intXY,2); %desired # int nodes (nIntNodes) is guaranteed to match with # intx after error checks
                    trial=trial+1;
                end
                
            case 's' %Stossell 2011 method
                
    
                %Form lookup table for initial filament between adjacent peripheral nodes
                
                %sort peripheral node angles
                [~, idx]=sort(periPol(1,:));
                
                %connect adjacent peripheral nodes
                adj_nodes1=idx;
                adj_nodes2=[idx(2:end) idx(1)];
                
                %Form lookup table for initial filaments from top half to
                %lower half (Cartesian second method)
                top_nodes=find(periXY(2,:)>er); %nodes on the top half
                bot_nodes=find(periXY(2,:) <-er); %code will need to account for when there are peripheral nodes at exactly the X-axis
                bot_nodes_match=NaN(1,length(top_nodes));
                
                %for each top node, find closest node between it and the
                %X-axis
                
                for ii=1:length(top_nodes)
                    
                    %% for each top node, find the the one bottom node that shares the same X-axis +/- small error
                    
                    matched_bot_nodes= find(abs(periXY(1,top_nodes(ii))-periXY(1,bot_nodes)) <= er);
                    
                    %if there are multiple matched bot nodes that share same
                    %X-axis as the top node, add a null match to be removed later.  Otherwise store the matched
                    %node
                    if length(matched_bot_nodes) > 1
                        bot_nodes_match(ii) = NaN;
                    else bot_nodes_match(ii) = matched_bot_nodes;
                    end
                end
                
                %remove top nodes with multiple matches, as well as
                %corresponding bot nodes as this means they are colinear 
                top_nodes(isnan(bot_nodes_match))=[];
                bot_nodes_match(isnan(bot_nodes_match)) = [];
                
                if nargin==0;	%display top/bottom filaments
                    filNodeLUT_periOnly_TB(1,:)=top_nodes;
                    filNodeLUT_periOnly_TB(2,:)=bot_nodes(bot_nodes_match);
                    figure; PrintNetwork_XY( periXY,filNodeLUT_periOnly_TB, 0, cellRadius*1.01, [], [], 'top/bottom nodes',true,nPeriMobile); %peripheral nodes only with peripheral filaments
                end
                
                %Form lookup table for initial filaments from left half to
                %right half (Cartesian second method)
                right_nodes=find(periXY(1,:)>er); %nodes on right half
                left_nodes=find(periXY(1,:)<-er);
                left_nodes_match=NaN(1,length(right_nodes));
                for ii=1:length(right_nodes)
                    
                    
                    %% for each right half node, find correpsonding left half node by reflecting across the Y-axis +/- small error
                    
                    matched_left_nodes = find( abs( periXY(2,right_nodes(ii)) - ( periXY(2,left_nodes ) ) ) <= er);
                    
                    %if there are multiple matched left nodes that share same
                    %Y-axis as the RIGHT node, add a null match to be removed later.  Otherwise store the matched
                    %node
                    if length(matched_left_nodes) > 1
                        left_nodes_match(ii) = NaN;
                    else left_nodes_match(ii) = matched_left_nodes;
                    end                  
                end
                
                %remove RIGHT nodes with multiple matches, as well as
                %corresponding LEFT nodes as this means they are colinear
                right_nodes(isnan(left_nodes_match))=[];
                left_nodes_match(isnan(left_nodes_match)) = [];
                
                if nargin==0	%display left/right filaments
                    filNodeLUT_periOnly_LR(1,:)=right_nodes;
                    filNodeLUT_periOnly_LR(2,:)=left_nodes(left_nodes_match);
                    figure; PrintNetwork_XY( periXY,filNodeLUT_periOnly_LR, 0, cellRadius*1.01, [], [], 'left/right nodes',true,nPeriMobile); %peripheral nodes only with peripheral filaments
                end
                
                %Combine top/bottom and left/right filaments into 1 initial
                %filaments LUT
                filNodeLUT_periOnly(1,:)=[adj_nodes1 top_nodes right_nodes];
                filNodeLUT_periOnly(2,:)=[adj_nodes2 bot_nodes(bot_nodes_match) left_nodes(left_nodes_match)];
                
                %Remove duplicate filaments
                filNodeLUT_periOnly_s=sort(filNodeLUT_periOnly);
                filNodeLUT_periOnly_s=union(filNodeLUT_periOnly_s',filNodeLUT_periOnly_s','rows'); %union of cols of matrix in LUT_new (i.e. removes duplicate cols). Sorted.
                filNodeLUT_periOnly=filNodeLUT_periOnly_s';
                
                if nargin==0    %display all filaments (no intersections left)
                    figure; PrintNetwork_XY( periXY,filNodeLUT_periOnly, 0, cellRadius*1.01, [], [], 'Peripheral nodes and peripheral filaments',true,nPeriMobile); %peripheral nodes only with peripheral filaments
                end
                
                if periNodePositioning == 's' && nPeriNodes == 12
                    disp ('ERROR GenerateNetwork: current code does not support finding intersections when there are 3 nodes per side (i.e., there is one node in top and one node in bottom hemispheres with the middle node on the X-axis directly)');
                    return
                end
                
                trial=0; nIntx=-1;
                while nIntx<nIntNodes %while the # int nodes not sufficient, loop
                    %Note: passing Error Checks guarantees that desired # int nodes will
                    %be found
                    
                    [intXY, intxFils]=FindIntersections_Schwarz(nIntNodes, periXY, filNodeLUT_periOnly);    %generate nIntNodes number of internal nodes at intersections of filaments
                    
                    nIntx=size(intXY,2); %desired # int nodes (nIntNodes) is guaranteed to match with # intx after error checks
                    trial=trial+1;
                end
        end
        
        if nargin==0
            disp(['# trials until sufficient intxs found: ',int2str(trial)]);
            figure; PrintNetwork_XY([intXY periXY],size(intXY,2)+filNodeLUT_periOnly, size(intXY,2), cellRadius*1.01, [], [], 'All nodes and peripheral filaments',true,nPeriMobile); %both peripheral nodes and internal nodes with peripheral filaments
        end
        
        %Add noise to intersection locations
        intXY=intXY+ intNodeNoise*cellRadius*rand(size(intXY) );
        
        [intPol(1,:), intPol(2,:)]=cart2pol(intXY(1,:),intXY(2,:));  %convert internal node positions to polar crds
        
        if ~isempty(intxFils) %if at least one intersection was found...
            filNodeLUT=UpdateFilaments(filNodeLUT_periOnly,intxFils,intXY,periXY,nFilaments);   %...normalize filNodeLUT and connect filaments to internal nodes
        else
            disp('ERROR (GenerateNetwork): there were no intersecting filaments found');
            filNodeLUT=filNodeLUT_periOnly;
        end
        
        if nargin==0
            disp(['# peripheral nodes: ',int2str(nPeriNodes)]);
            disp(['initial # peri filaments: ',int2str(nPeriFilaments)]);
            disp(['# internal nodes: ',int2str(nIntNodes)]);
            disp(['desired # filaments: ',int2str(nFilaments)]);
            disp(['actual # filaments: ',int2str(size(filNodeLUT,2))]);
        end
end

%plot network after normalizing/updating filNodeLUT
if nargin==0
    %     figure; PrintNetwork_XY([intXY periXY],filNodeLUT,nIntx,cellRadius*1.01,[],[],'Final network',true);
    nodeXY=[intXY periXY];
    figure;     polar(1,cellRadius*1.01);
    hold on; PlotLabeledNetwork(filNodeLUT,nodeXY,nodeXY(:,1:nIntNodes),nodeXY(:,nIntNodes+1:end),nIntNodes,true,[],nPeriMobile);
    title('Final network w/ noise (optional)');
end
end

%APPENDIX
% function []=Appendix()
%
% %Polar method
% bot_nodes_p=NaN(1,length(top_nodes));
% for ii=1:length(top_nodes)
%     bot_nodes_p(ii)=find( abs( periPol(1,:) - ( 2*pi-periPol (1,top_nodes(ii) ) ) ) <= er ); %for each top half node, find correpsonding bottom half node by reflecting across the X-axis +/- small error
% end
% filNodeLUT_periOnly_p(1,:)=top_nodes;
% filNodeLUT_periOnly_p(2,:)=bot_nodes_p;
% PrintNetwork_XY( periXY,filNodeLUT_periOnly_p, 0, cellRadius, [], [], 'polar method',true); %peripheral nodes only with peripheral filaments
%
% %Cartesian method
% filNodeLUT_periOnly_1(1,:)=top_nodes; %origin nodes are the top nodes
% for ii=1:length(top_nodes)                 %fill in destination nodes as bot nodes
%     bn_sameX= abs(periXY(1,top_nodes(ii))-periXY(1,bot_nodes)) < er ; %for each top node, find the bottom nodes that shares the same X-axis
%     filNodeLUT_periOnly_1(2,ii)=bot_nodes(bn_sameX);
% end
% figure; PrintNetwork_XY( periXY,filNodeLUT_periOnly_1, 0, cellRadius, [], [], 'cartesian method',true); %peripheral nodes only with peripheral filaments
%
% %plot network before normalizing/updating filNodeLUT
% if nargin==0
%     PrintNetwork_Polar([],0,1,1,[],periPol,filNodeLUT_periOnly,[],100,1);   %only plots peripheral nodes and their connecting filaments
%     if ~isempty(intPol); polar(intPol(1,:),intPol(2,:),'ko'); end  	%plots internal nodes as circles
%     title('filNodeLUT_periOnly','Interpreter','none');
%     for ii=1:size(filNodeLUT_periOnly,2)
%         text(periXY(1,filNodeLUT_periOnly(:,ii))+10*randi([-1,1]),periXY(2,filNodeLUT_periOnly(:,ii))+10*randi([-1,1]),int2str(ii)); %label both ends of initial filaments (between only peri nodes)
%     end
% end
%
% %plot network after normalizing/updating filNodeLUT
% if nargin==0
%     [allXY(1,:), allXY(2,:)]=pol2cart([intPol(1,:) periPol(1,:)],[intPol(2,:) periPol(2,:)]);
%     PrintNetwork_Polar([],0,1,1,intPol,periPol,filNodeLUT,[],100,1);
%     title('filNodeLUT');
%     for ii=1:size(filNodeLUT,2)
%         text(allXY(1,filNodeLUT(:,ii))+10*randi([-1,1]),allXY(2,filNodeLUT(:,ii))+10*randi([-1,1]),int2str(ii));    %label both ends of final filaments (between peri and int nodes)
%     end
%     end
% end