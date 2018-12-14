%Is there redundancy between this and FindPairedElement function for
%finding the columns?

function [nodeFilLUT,nFilsPerNode,nFilsPerOrgn] = PopulateNodeFilLUT(nSearched,filNodeLUT)
%PopulateNodeFilLUT: populate LUTs for both which and how many filaments
%are connected to each of the first "nSearched" nodes in filNodeLUT. 

%       OUPUTS:
%           -nodeFilLUT: nSearched-by-nSearched matrix describing "node->connected filaments" LUT to look up connected filaments
%               -rows describe fils, cols describe nodes
%           -nFilsPerNode: how many connected filaments total per node

%           -nFilsPerOrgn: how many filaments with node as origin node
%           ("top" node vs. dest node/"bottom" node)

%		NOTE: both arrays/LUT contain mobile nodes only (i.e., not fixed peripheral nodes)

%		NOTE: code will search through FIRST "nSearched" nodes in
%		"filNodeLUT"...note that the mobile nodes come first, then
%		non-mobile nodes (i.e., to search through non-mobile nodes as well,
%		simply increase nSearched)

%12/01/18: added test cases

if nargin==0
   sample=1;
   switch sample
       case 1
           %3 filaments
           nSearched=2;
           filNodeLUT=[1 3; 1 2; 4 5]';
   end
end

%Col rep nodes; for each col, rows rep filaments connected to node ("column number matched node number")
%May want to change nSearched. (?)
nodeFilLUT = zeros(nSearched);  %nSearched-by-nSearched zero matrix

%Row-indexed nodes. Values rep # filaments connected to node (rep length of columns in nodeFilLUT)
nFilsPerNode = zeros(nSearched,1);

%Default [nAllNodes,1]?? Row-indexed nodes. Values rep # filaments that
%have node as "origin" node (vs. "destination" node)
nFilsPerOrgn = zeros(nSearched,1);

for ii = 1:nSearched	%only searches filNodeLUT for first "nSearched" nodes
    origFils = find(filNodeLUT(1,:)==ii);    %vector of filaments containing "origin" node ii
    destFils = find(filNodeLUT(2,:)==ii);    %vector of filaments containing destination" node ii
    %     if ~isempty(origFils) || ~isempty(destfils)         %If Node ii is contained by any filaments...
    %         nodeFilLUT(:,ii)=cat(1,origFils',destFils');	%fill in filaments with ii as "origin" first then fill in filaments with Node ii as "destination"
    %         nFilsPerOrgn(ii,1)=length(origFils);       %keeps track of how many filaments node ii is the origin in
    %         nFilsPerNode(ii,1)=length(origFils)+length(destFils); %vector storing how many filaments node ii is in (both origin and dest)
    %     end
    nOrgn = length(origFils);    %# of times node ii appears as the "origin" node
    nDest = length(destFils);    %# of times node ii appears as the "destination" node
    totFils = nOrgn+nDest;		 %total number of filaments that contain node ii
    if totFils > 0               %If Node ii is contained by any filaments...
        
        %In column ii, fill in filaments with Node ii as "origin" first, then fill in filaments with Node ii as "destination"
        nodeFilLUT(1:totFils,ii)=cat(1,origFils',destFils');
        
        nFilsPerOrgn(ii,1) = nOrgn;    %keeps track of how many filaments node ii is the origin in
        nFilsPerNode(ii,1) = totFils;    %vector storing how many filaments node ii is in (both origin and dest)
    end
end
end

