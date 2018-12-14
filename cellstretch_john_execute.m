%11/09/09
%This file contains functions:
%1. cellstretch_john_execute
%2. BreakFilamentsByStrain
%3. StretchNodes
%4. CopyCodeToDir (not implemented)

function [lastCycle]=cellstretch_john_execute(SAVE_DIR,whichRun,numCycles,datasetDescription,doElasticNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,stretchMag,breakMethod,breakFactor,doRotate,cellRadius,intNodeNoise,elastMod,nIntNodes,nPeriNodes,nFilaments,delX,relaxThresh,iMax,doSaveNetworkBeforeBreak,figOptions)

%12/01/23: changed from tracking forces (which depends on elastMod) to
%tracking strains.  Subsequently changed relaxation threshold from a sum of
%squared forces to a sum of squared displacements [PENDING]

%12/01/17: changed the way ministretch is performed from 1:delX:stretchMag
%to delX:delX:stretchMag.  Track the amt of stretch (currStretch) rather
%than the number of ministretches (miniStretchCount)

% disp(['[Run ' int2str(whichRun) '] Begin run at ' GetCurrentTime() ]);
disp(['[R' int2str(whichRun) '] Folder name: ' SAVE_DIR]);
disp(['[R' int2str(whichRun) '] Dataset description: ' datasetDescription])
disp(['[R' int2str(whichRun) '] delX = ' num2str(delX) ', relaxThresh = ' num2str(relaxThresh)]);
rand('state',sum(100*clock));   %resets random generator to a different state each time

%GENERATE INITIAL NETWORK
disp(['[R' int2str(whichRun) '] ' GetCurrentTime() ' generating network...' ]);
[intPol_0,periPol_0,nPeriMobile,filNodeLUT]=GenerateNetwork(nIntNodes,nPeriNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,nFilaments,cellRadius,intNodeNoise,doRotate);
disp(['[R' int2str(whichRun) '] ' GetCurrentTime()  ' finished generating network' ]);

nodePol_0=[intPol_0 periPol_0];   %Left-to-Right: internal nodes first, then peripheral nodes.  "_0" denotes initial position.

nMobile=nIntNodes+nPeriMobile;  %mobile # nodes

nNodes = nPeriNodes + nIntNodes; %total # nodes, =sizeall

[nodeXY(1,:),nodeXY(2,:)]=pol2cart(nodePol_0(1,:),nodePol_0(2,:)); %nodeXY: mobile nodes first, then fixed nodes

%PRINT/SAVE LABELED NETWORK
if figOptions(1); figure; else figure('Visible','off'); end
TITLE=['Run',int2str(whichRun),'_Cycle0_LabeledNetwork'];
polar(pi,cellRadius+delX);
PlotLabeledNetwork(filNodeLUT,nodeXY,nodeXY(:,1:nIntNodes),nodeXY(:,nIntNodes+1:end),nIntNodes,true,TITLE);
if figOptions(2); print('-dtiff',[SAVE_DIR TITLE]); end

%PRINT/SAVE INITIAL NETWORK PRIOR TO STRETCHING
if figOptions(1); figure; else figure('Visible','off'); end
TITLE = ['Run',int2str(whichRun),'_Cycle0_InitialNodes'];
PrintNetwork_XY(nodeXY,[]        ,nIntNodes,cellRadius+delX,[],[],TITLE,false,nPeriMobile); %FIGURE 1 - focal adhesions/peripheral nodes ONLY (no filaments)
if figOptions(2); print('-dtiff',[SAVE_DIR TITLE]); end

if figOptions(1); figure; else figure('Visible','off'); end
TITLE =(['Run',int2str(whichRun),'_Cycle0_InitialNetwork']);
PrintNetwork_XY(nodeXY,filNodeLUT,nIntNodes,cellRadius+delX,[],[],TITLE,false,nPeriMobile);	 %FIGURE 2 - entire network including both nodes & filaments
if figOptions(2); print('-dtiff',[SAVE_DIR TITLE]); end

%SAVE VARIABLES BEFORE STRETCHING
save([SAVE_DIR TITLE]);
nodeXY_0=nodeXY;    %saves original positions of nodes

%===***MAIN: stretches filaments by delX for numCycles cycles, relax to equilibrium, find which strains are > threshold
breakCounter = zeros(1,numCycles); %tracks # filaments broken on each cycle
for whichCycle = 1:numCycles
    if doElasticNodes %Kathy's original method that restores the initial positions of the interior nodes every cycle (i.e. elastic nodes)
        
        %Restore nodeXY using saved nodeXY;
        nodeXY=nodeXY_0; %RESETS the locations of ALL nodes
        
    else %only restores fixed nodes but does not reset positions of mobile nodes
        
        [nodeXY(1,nMobile+1:nNodes), nodeXY(2,nMobile+1:nNodes)] = pol2cart(nodePol_0(1,nMobile+1:nNodes),nodePol_0(2,nMobile+1:nNodes)); %only restore "fixed" perimeter nodes
        
    end
    
    %For mobile nodes only: populate node->which filaments LUT and populate node->how many connected filaments
    [connectivity.nodeFilLUT,connectivity.nFilsPerNode,connectivity.nFilsPerOrgn] = PopulateNodeFilLUT(nMobile,filNodeLUT);
    %"PopulateNodeFilLUT" outputs is only used for SumStrains function. SumStrains is called every "mini-cycle"
    
    [oldFilLen,~]=FindFilLengthsAndAngles(nodeXY,filNodeLUT); %Find the original lengths of the filaments (e.g. after previous cycle's relaxing) before any stretching
    
    %==********************************************************************
    %==Stretch fixed nodes and relax mobile nodes
    %==********************************************************************
    %     miniStretchCounter = 0;             %counter for # "micro-stretches".
    %      Don't need to "count" anymore since delX can be any number and is
    %      more useful to track total stretch instead of number of mini
    %      stretches.
    
    %==Start mini-stretches & relax
    for currStretch = delX:delX:stretchMag    %"for loop" iterates (stretchMag-delX)/delX times (9 (10??) by default). Total stretch is stretchMag
        
        %STRETCH FIXED NODES
        %-color filaments such that the negative strains can be
        %visualized??
        
        disp(['[R' int2str(whichRun) ' C' int2str(whichCycle) ' S' num2str(currStretch)  '%] ' GetCurrentTime() ' stretching....'])
        
        switch stretchMethod
            case 'p' %peripheral node stretch method (default)
                %Uniaxially stretchs the chosen nodes outwards in X-axis
                %direction by delX amount (<---->)
                stretchNodes=nMobile+1:nNodes;
                nodeXY(1,stretchNodes) = nodeXY(1,stretchNodes).*(1+delX/100); %row 1 of nodeXY is X-axis
                
            case 'i' %internal node stretch method (use this for shearing in the future)
                %stretchNodes=internal nodes
                
            case 's' %Stossel method
                %Uniaxially stretchs the chosen nodes leftwards in X-axis
                %direction by delX amount (-->)
                
%                 stretchNodes=nIntNodes+FindDataSubset(periPol_0,pi/4,3/4*pi,1,'column data'); %upper_quadrant_nodes;
                stretchNodes=nIntNodes+FindDataSubset(periPol_0,(90-vertOffset)*pi/180,(90+vertOffset)*pi/180,1,'column data'); %upper_quadrant_nodes;
                nodeXY(1,stretchNodes)=nodeXY(1,stretchNodes)+delX/100 *cellRadius *2;  %"*2" factor because Stossel stretch is calclated as dx/h where dx=actual stretch, h=diameter 
            otherwise
                disp(['ERROR (cellstretch_john_execute): stretchMethod ' stretchMethod ' not recognized']);
                return;
        end
        
        %         miniStretchCounter = miniStretchCounter+1;
        
        %PRINTS NETWORK AFTER INCREMENTAL STRETCH
        if mod(currStretch,1)==0  %only output figures when stretch % is an integer
            if figOptions(1); figure; else figure('Visible','off'); end
            TITLE=(['Run' int2str(whichRun) '_Cycle' int2str(whichCycle) '_StretchPct' num2str(currStretch) '_afterStretch.tif']);
            PrintNetwork_XY(nodeXY,filNodeLUT,nIntNodes,cellRadius+delX,[],[],TITLE,false,nPeriMobile);
            if figOptions(2); print('-dtiff',[SAVE_DIR TITLE]); end
        end
        
        
        %RELAX MOBILE NODES IN NETWORK
        %calculate forces of each node and adjust location until all nodes are in equilibrium (summation of forces is ~zero)
        %GS_thresh=0.1;
        disp(['[R' int2str(whichRun) ' C' int2str(whichCycle) ' S' num2str(currStretch)  '%] ' GetCurrentTime() ' start relaxing....'])
        
        %-New function: [nodeXY_]=RelaxNetwork(nodeXY,filNodeLUT,oldFilLen,GS_thresh,connectivity,nFirst=nMobile,);
        
        %         fSquare = intmax;                        %represents summed squared forces (strain*E)^2 of all filaments
        %         fSquareAll = zeros(iMax,1);         %stored fSquare values for each iteration
        sSquare= intmax;    %represents summed squared strains (sum(strain^2)) of all filaments
        sSquareAll = zeros(iMax,1);  %stored sSquare values for each iteration
        
        ii=1;                    		%iteration counter, must be < iMax
        strainFactor = 1;             	%strainFactor is the factor that strain is multiplied by to convert it to step size
        
        %         while (fSquare >= relaxThresh*elastMod)                %ii++ every loop. Escape condition if ii > iMax
        while (sSquare >= relaxThresh)                %ii++ every loop. Escape condition if ii > iMax
            
            %11/08/02
            [strcFilLen,strcFilAng]=FindFilLengthsAndAngles(nodeXY,filNodeLUT);   %find new filament lengths and angles
            
            %"Strain" is a vector of magnitudes for all filaments
            Strain = double((strcFilLen-oldFilLen)./oldFilLen);     %strcFilLen contains the "new" filament lengths. oldFilLen is the lengths before any stretch this cycle.
            
            %sum strains on the first nMobile nodes (the mobile nodes).
            %Note: nodeStrainXY contains only strains for only mobile nodes
            nodeStrainXY = SumStrains(Strain,strcFilAng,connectivity);
            
            %adjust node positions (nodeXY) based on Strains
            %             [nodeXY,fSquare]=AdjustMovableNodes(nMobile,nodeXY,strainFactor,nodeStrainXY,elastMod);
            %             [nodeXY,fSquare]=AdjustMovableNodes(nodeXY(:,1:nMobile),strainFactor,nodeStrainXY,elastMod);
            nodeXY(:,1:nMobile)= nodeXY(:,1:nMobile) + strainFactor.*nodeStrainXY;   %adjust node positions in the horizontal direction
            
            
            %changed sumForceSquare to sumStrainSquare and removed elastMod from the input parameters since it can easily
            %be multiplied by the sumStrainSquare to become sumForceSquare outside the
            %function.
            sSquare = sum(sqrt(nodeStrainXY(1,:).^2 + nodeStrainXY(2,:).^2));   %sum squared strains for the mobile nodes in current iteration
            
%             fSquare=sumStrainSquare*elastMod;
%             fSquareAll(ii)=fSquare;
            
            sSquareAll(ii)=sSquare; %stores sum squared strain for the current iteration of relaxation
            
            if ii >= iMax
                break
            elseif ii > 1
%                 strainFactor=AdjustStrainFactor(strainFactor,fSquareAll(ii-1),fSquareAll(ii));
                strainFactor=AdjustStrainFactor(strainFactor,sSquareAll(ii-1),sSquareAll(ii));
            end
            ii = ii+1;  %iteration counter
            
            %OUPUTS DURING RELAXATION OF MOBILE NODES
            if mod(ii,iMax*0.05)==0 %print network when at certain iterations
                STATUS=['[R' int2str(whichRun) ' C' int2str(whichCycle) ' e' num2str(currStretch) '%] ' GetCurrentTime() ' ii=' int2str(ii) ' max abs strain=' num2str(max(max(abs(nodeStrainXY)))) ' sSquare=' num2str(sSquare)];
                disp(STATUS);
                
                %PRINT NETWORK WITH STRAINS OF FILAMENTS AND NODES
                if figOptions(1); figure; else figure('Visible','off'); end
                
                TITLE=['Run' int2str(whichRun) '_Cycle' int2str(whichCycle) '_e' num2str(currStretch) '_ii' int2str(ii) '_maxAbsStrain' num2str(max(max(abs(nodeStrainXY)))) '_sSquare' num2str(sSquare) '.tif'];
                PrintNetwork_XY(nodeXY, filNodeLUT, nIntNodes, cellRadius+delX, [], [], TITLE, [false false false true true], nPeriMobile, nodeStrainXY, Strain);
                
                if figOptions(2) && whichCycle==1 && whichRun==1; print('-dtiff',[SAVE_DIR TITLE]); end
                
                %PRINT STRAINS COLUMN-INDEXED BY NODES
                if figOptions(1); figure; else figure('Visible','off'); end
            end
        end
        %/End RELAXATION OF MOBILE NODES
        
        disp(['[R' int2str(whichRun) ' C' int2str(whichCycle) ' S' num2str(currStretch)  '%] ' GetCurrentTime() ' finished relaxing'])
        
        %PRINTS NETWORK AFTER RELAXATION OF MOBILE NODES
        if mod(currStretch,1)==0  %only output figures when stretch % is an integer
            if figOptions(1); figure; else figure('Visible','off'); end
            TITLE=(['Run' int2str(whichRun) '_Cycle',int2str(whichCycle) '_StretchPct' num2str(currStretch) '_postRelax_ii' int2str(ii) '.tif' ]);
            PrintNetwork_XY(nodeXY,filNodeLUT,nIntNodes,cellRadius+delX,[],[],TITLE,false,nPeriMobile);
            if figOptions(2) && whichCycle==1 && whichRun==1; print('-dtiff',[SAVE_DIR TITLE]); end
        end
        
        disp(['[R' int2str(whichRun) ' C' int2str(whichCycle) ' S' num2str(currStretch)  '%] ' GetCurrentTime() ' term iter: ' int2str(ii) ', strain factor ' num2str(strainFactor)]);
        if ii == iMax	%if has not relax under iMax cycles, ...
            break	%...then abandon minicycle stretching
        end
        
        %Save network after relaxation
        name =['Run' int2str(whichRun) '_Cycle' int2str(whichCycle) '_StretchPct' num2str(currStretch) '_postRelax.mat'];
        save([SAVE_DIR name]);
        
    end %==End Mini-stretch & relax
    %==********************************************************************
    %==/End stretch fixed & relax mobile nodes
    %==********************************************************************
    
    %12/01/18: is this necessary when the last "postRelax" network is the same?
%     if doSaveNetworkBeforeBreak==1
%         name =['Run' int2str(whichRun) '_Cycle' int2str(whichCycle) '_StretchPct' num2str(currStretch) '_afterTotalStretch_beforeBreak']; 
%         save([SAVE_DIR name]);  %SAVES VARIABLES BEFORE BREAKAGE
%     end
    
    if ii == iMax	%if iMax hit, exit outer cycle loop and return to main function to attempt next filament network
        break
    end
    
    %==********************************************************************
    %==Determine which filaments to break and break them
    %==********************************************************************
    
    disp(['[R' int2str(whichRun) ' C' int2str(whichCycle) '] ' GetCurrentTime() ' breaking and replacing...'])
    
    switch breakMethod
        case 'k' %Kathy
            %[filNodeLUT_afterReplacement,fBreak] = BreakAndReplaceFilaments_Kathy(Strain,breakFactor,filNodeLUT);
            [filNodeLUT_break,fBreak]=BreakFilamentsByStrain(filNodeLUT,Strain,breakFactor);
            filNodeLUT_afterReplacement=GenerateRandomFilaments(length(fBreak),nNodes,filNodeLUT_break);
            filGenCount=size(fBreak,2);
            
        case 'c' %Critical concentration
            %Critical concentration=sum of the lengths of all filaments in
            %the fully stretched network
            conc_crit=sum(FindFilLengthsAndAngles(nodeXY,filNodeLUT))*breakFactor;
            
            disp(['R[' int2str(whichRun) ' C' int2str(whichCycle) '] critical conc: ',num2str(conc_crit)]);
            
            %Proposed network after breakage/reform
            [filNodeLUT_break,fBreak]=BreakFilamentsByStrain(filNodeLUT,Strain,breakFactor);
            currentSum=-1;
            
            %Iteratively reform filaments
            filGenCount=0;
            while currentSum<conc_crit
                %should new filaments be generated across peripheral nodes
                %AND internal nodes or only peripheral nodes?
                filNodeLUT_break=GenerateRandomFilaments(1,nNodes,filNodeLUT_break);
                currentSum=sum(FindFilLengthsAndAngles(nodeXY,filNodeLUT_break));
                filGenCount=filGenCount+1;
            end
            disp(['R',int2str(whichRun),' C',int2str(whichCycle),' filaments broken: ',int2str(size(fBreak,2)),', filaments generated: ' int2str(filGenCount)]);
            
            filNodeLUT_afterReplacement=filNodeLUT_break;
            
        otherwise
            disp(['ERROR (cellstretch_john_execute): break method ' breakMethod ' not recognized'])
    end
    
    disp(['[R' int2str(whichRun) ' C' int2str(whichCycle) '] ' GetCurrentTime() ' finished breaking and replacing'])
    
    %DETERMINE FILAMENTS FOR BREAKAGE
    %     [fBreak]=FindBreakFilaments(breakMethod,breakFactor,Strain); %check which filaments break
    
    %REPLACE CHOSEN FILAMENTS WITH NEW ONES
    %     filNodeLUT_afterReplacement=ReplaceFilaments(intNodePositioning,fBreak,filNodeLUT,nodeXY,nIntNodes); %create new filaments, connected to different nodes, to replace these. Other LUTs will be updated next cycle
    
    %PRINTS NETWORK AFTER BREAKAGE
    breakCounter(whichCycle)=length(fBreak);
    if figOptions(1); figure; else figure('Visible','off'); end
    TITLE=(['Run' int2str(whichRun) '_Cycle' int2str(whichCycle) '_terminal']);
    PrintNetwork_XY(nodeXY,filNodeLUT_afterReplacement,nIntNodes,cellRadius+delX,[],[],[TITLE '_w_',int2str(size(fBreak,2)),'_fils_broken_',int2str(filGenCount),'_fils_generated'],false,nPeriMobile);
    doPrint=1:4:numCycles; %if whichCycles is one of the cycles in doPrint, then print
    if figOptions(2) && sum(whichCycle==doPrint); print('-dtiff',[SAVE_DIR TITLE]); end
    
    %Save completed network
    name =['Run',int2str(whichRun),'_Cycle',int2str(whichCycle),'_complete'];
    save([SAVE_DIR name]);
    
    filNodeLUT=filNodeLUT_afterReplacement;
    
end %/end outer cycle loop

lastCycle=whichCycle;
save([SAVE_DIR 'last_cycle'],'lastCycle'); %saves only the last_cycle variable
clear variables;
end

function [filNodeLUT_break,fBreak]=BreakFilamentsByStrain(filNodeLUT,Strain,breakFactor)
%BreakFilamentsByStrain: breaks filaments probabilistically such that
%filaments with heaviest strain are more likely to break
%-Note: if breakFactor is larger (i.e., >1), then it is less likely
%filaments will break

nFilaments=size(Strain,2);
normStrain = abs(Strain)./max(abs(Strain));
fBreak = find(normStrain>breakFactor*rand(1,nFilaments)); %determine which filaments break and store indices in filamentsToBreak

filNodeLUT_break=filNodeLUT;
filNodeLUT_break(:,fBreak)=[];
end


function []=CopyCodeToDir(SAVE_DIR)
%this function must copy used code files to the current data directory

end