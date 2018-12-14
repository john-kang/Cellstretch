function [filamentsN]=UpdateFilaments(filaments,intxFils,intXY,periXY,nMax)
%UpdateFilaments: given the LUT for filaments->nodes ("filaments"); how the
%filaments intersect (intxFils); the XY-crds for the peripheral (periXY)
%and internal nodes (intXY); and what the final filament count (nMax)
%should be, this func generates a new LUT (filamentsN) that's normalized to
%account for new interior nodes.  The number of filaments in filamentsN
%will be as close to nMax as possible given the original network
%intersections but will not be greater than nMax.

%11/07/23: major revision confirmed to work properly.  When nMax (max
%number of total filaments) is altered, the filamentsN LUT is altered
%correctly.  Confirmed with plotting using function "PlotLabeledNetwork".

%11/07/19: major revision should be working; proper filament number is
%generated.  Revising code to improve reliability.

%11/07/16: major revision is largely working.  Need to implement nMax in
%order to generate the proper number of filaments

%11/07/12: MAJOR REVISION -code did not work properly as intended...
%-removed nIntNodes as input, now taken directly from intxFils

%11/07/08: normalizes by adding nIntNodes to "filaments".  This accounts
%for the fact that "filaments" initially only assumes peripheral nodes
%exist but "filamentsN" assumes that internal nodes/intersections exist as
%well and thus needs the offset

%11/06/29: initial version

if nargin==0
    clear all; close all;
    disp('Start UpdateFilaments');
    example=3;
    switch example
        case (1) %5 peri nodes, 2 intersections
            filaments=[1,2,2;3,4,5]; %3 fils
            intxFils=[1,1;3,2]; %2 intersections
            periXY=[80.9016994374947,-30.9016994374948,-100,-30.9016994374948,80.9016994374947;58.7785252292473,95.1056516295154,1.22464679914735e-14,-95.1056516295154,-58.7785252292473;];
            intXY=[11.8033988749895,-30.9016994374948;36.3271264002681,22.4513988289793;];
            nMax=4;
        case (2) %
            filaments=[1,1,1,1,3,3,4,5,5,8;6,7,8,9,5,6,9,6,10,10;];
            intxFils=[3,4,1;6,10,5;];
            periXY=[6.12323399573677e-15,-1.83697019872103e-14,95.1056516295154,58.7785252292473,-58.7785252292473,-95.1056516295154,-95.1056516295154,-58.7785252292474,58.7785252292473,95.1056516295154;100,-100,30.9016994374947,80.9016994374947,80.9016994374947,30.9016994374948,-30.9016994374947,-80.9016994374947,-80.9016994374948,-30.9016994374948;];
            intXY=[-22.4513988289793,47.5528258147577,-36.3271264002681;30.9016994374948,-46.3525491562421,73.6067977499790;];
            nMax=20;
        case (3) %
            filaments=[1,2,4,5,5,6,6,7,8,9;7,7,9,6,8,7,9,10,10,10;];
            intxFils=[1,2,5,3;7,5,8,9;];
            periXY=[6.12323399573677e-15,-1.83697019872103e-14,95.1056516295154,58.7785252292473,-58.7785252292473,-95.1056516295154,-95.1056516295154,-58.7785252292474,58.7785252292473,95.1056516295154;100,-100,30.9016994374947,80.9016994374947,80.9016994374947,30.9016994374948,-30.9016994374947,-80.9016994374947,-80.9016994374948,-30.9016994374948;];
            intXY=[-65.7163890148917,-58.7785252292474,-58.7785252292474,58.7785252292473;9.54915028125265,-57.2949016875157,-30.9016994374947,-42.7050983124842;];
            nMax=20;
    end
    PlotLabeledNetwork(filaments,periXY,intXY,periXY,0);     %Plots initial network
end

%ERROR CHECKS
if size(intxFils,2)~=size(intXY,2)
    disp('ERROR (UpdateFilaments): Number of intersections and intersection coordinates not consistent'); end
if max(max(filaments))>size(periXY,2)
    disp('ERROR (UpdateFilaments): filament definitions not consistent with number of peripheral points'); end
if size(intXY,2)>nchoosek(size(periXY,2),2)
    disp('ERROR (UpdateFilaments): too many intersections for given number of peripheral points'); end
if max(max(intxFils))>size(filaments,2)
    disp(['ERROR (UpdateFilaments): number of filaments ',int2str(size(filaments,2)),' not consistent with how filaments intersect (max fil ID is ',int2str(max(max(intxFils))),')']); end
if size(filaments,2)>=nMax
    disp(['ERROR (UpdateFilaments): number of initial filaments (',int2str(size(filaments,2)),') is greater than or equal to max number of filaments (',int2str(nMax),').  No new filaments will be created']); end

filamentsN=filaments+size(intXY,2); %offset the peripheral nodes because the internal nodes have lower indices

%For each intersection, find two intersecting filaments
for intx=1:size(intXY,2)
    
    % For each filament of the intx
    for which=1:2
        thisFil=intxFils(which,intx); %"this filament" is the currently analyzed filament in the smallest subloop
        if nargin==0; disp(['...checking fil ',int2str(thisFil)]); end
        
        % If this filament has not been analyzed yet (vertices are not NaN)
        % and there is room for more filaments (nLeft>0)
        nFilsRemove=sum(mean(isnan(filamentsN))); %the number of original filaments that were analyzed for intersections.  These filaments are removed at the end and thus need to be subtracted from the total number of filaments to calculate the net number of filaments.
        nFils=size(filamentsN,2)-nFilsRemove; %the net number of filaments (taking into account filaments that will be removed at the end)
        nLeft=nMax-nFils; %the net number of filaments left to be created
        if ~isnan(filamentsN(:,thisFil)) & nLeft > 0
            
            %Find the intersections that intersect this filament
            [~,intxID]=FindPairedElements(intxFils,thisFil); %note that length(intxID)+1 is the num of fils that the original fil is broken into, but since the original filament is "lost", length(intxID) is the net num of fils generated
            %Note that intxID is essentially the ID of the interior nodes
            
            %Find Euclidean distance between intersection and one end of
            %this fil
            origin=filaments(1,thisFil);   %one end of this fil is treated as "origin"
            originXY=periXY(:,origin);   %XY crds for one end of this peri fil
            
            %if the net # fils generated (length(intxID)) is < the # fils
            %remaining, then store all intersection point
            if length(intxID) <= nLeft
                intxXYs=intXY(:,intxID);   %XY crds for all points that intx this fil
                
                %else (if the net # fils generated (intxID) exceeds or
                %equals the max # of new filaments), then truncate (i.e.,
                %remove one or more intx) and stop adding new filaments
            else
                intxXYs=intXY(:,intxID(1: nLeft));    %truncated XY crds for all points that intx this fil
                intxID=intxID(1:nLeft);
            end
            
            intx_dist=pdist([originXY intxXYs]'); %find dist b/w origin and other intx to this fil.
            %For m-by-n matrix, the element ordering is: [(2,1) (3,1) ,...,
            %(m,1),(3,2),...,(m,2),...,(m,m-1)).
            intx_dist=intx_dist(1:size(intxXYs,2));  %only interested in the first m pairings
            
            %rank distances
            [~,idx]=sort(intx_dist); %i.e., [[1 2 3],[2 3 1]]=sort([3 1 2]);
            
            %using rankings and location of each intersection, create an
            %array of new filaments
            add=zeros(2,length(intxID)+1); %# new filaments is # intx + 1
            add(:,1)=[filamentsN(1,thisFil);intxID(idx(1))]; %first new filament
            add(:,end)=[intxID(idx(end)),filamentsN(2,thisFil)]; %last new filament
            for ii=2:length(intxID)
                add(:,ii)= [intxID(idx(ii-1)) intxID(idx(ii)) ]'; %fill in the intermediate filaments
            end
            
            %concactenate additional filaments to end of norm filaments
            filamentsN=[filamentsN add];
            
            %indicate that this filament has been checked now to prevent
            %analyzing again
            filamentsN(:,thisFil)=NaN;   %make it so this filament cannot be "accessed" again. Should be a redudant mechanism to prevent repeat access along with the "checked" LUT.
            
        end
    end
end

%remove all the original filaments (vertices replaced with NaN in the LUT)
%that had been segmented into smaller filaments.  Then sort column-wise.
[~,colsToRemove]=find(isnan(filamentsN));
filamentsN(:,colsToRemove)=[];
filamentsN=sort(filamentsN);

if nargin==0
    figure;
    PlotLabeledNetwork(filamentsN,[intXY periXY],intXY,periXY,size(intXY,2));
end

if size(filamentsN,2)~=nMax
    disp(['Warning (UpdateFilaments): final number of filaments (',int2str(size(filamentsN,2)),') is not equal to desired number of filaments (',int2str(nMax),')']);
end
end
