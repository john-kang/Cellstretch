function [nodeStrainXY] = SumStrains(Strain, strcFilAng, connectivity)

%SumStrains: Adds strains on nodes (all internal nodes and certain peripheral nodes)
%Sum the strains on all mobile nodes (internal node and mobile peripheral node).  Strain proportional force by Hooke's Law
%Part of the Gauss Seidel algorithma

%nodeStrainXY sign "points" from origin node to destination node

nodeFilLUT=connectivity.nodeFilLUT; %for each mobile node, what filaments are connected
nFilsPerNode=connectivity.nFilsPerNode; %for each mobile node, how many filaments are connected
nFilsPerOrgn=connectivity.nFilsPerOrgn; %for each mobile node, how many filaments are connected to it at the "origin" as opposed to destination (must be <nFilsPerNode)
nMobile=length(nFilsPerNode);   %the number of mobile nodes is equal to the length of parameter nFilsPerNode, which is row-indexed mobile nodes

%ERROR CHECKS
if length(nFilsPerNode)~=length(nFilsPerOrgn)
    disp(['ERROR (SumStrains): nFilsPerNode (length ' int2str(nFilsPerNode) ') and nFilsPerOrgn length(' length(nFilsPerOrgn) ') not consistent']);
end

nodeStrainXY = NaN(2,nMobile);
for jj = 1:nMobile                 %mobile nodes listed first.  Excludes "stationary" peripheral nodes
    n1 = nFilsPerOrgn(jj,1);       %# filaments with node jj as origin
    n12 = n1+1;                    %divider between origin node and rest of nodes
    n2 = nFilsPerNode(jj,1);       %total # filaments w/ node jj
    
    %Sum strain of filaments with node jj as "origin" first, then...
    nodeStrainXY(1,jj) = sum(Strain(nodeFilLUT(1:n1,jj)).*cos(strcFilAng(1,nodeFilLUT(1:n1,jj))));
    nodeStrainXY(2,jj) = sum(Strain(nodeFilLUT(1:n1,jj)).*sin(strcFilAng(1,nodeFilLUT(1:n1,jj))));
    
    %...sum total strain of filaments with node jj as "destination".  
    %NOTE: This part is separate due to needing to add "pi"
    nodeStrainXY(1,jj) = nodeStrainXY(1,jj) + sum(Strain(nodeFilLUT(n12:n2,jj)).*cos(strcFilAng(1,nodeFilLUT(n12:n2,jj))+pi)); %"pi" makes up for offset due to axis origin at "origin" node
    nodeStrainXY(2,jj) = nodeStrainXY(2,jj) + sum(Strain(nodeFilLUT(n12:n2,jj)).*sin(strcFilAng(1,nodeFilLUT(n12:n2,jj))+pi));
end

end

