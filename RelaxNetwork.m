function [nodeXY_relaxed,msg] = RelaxNetwork(nodeXY,filNodeLUT,oldFilLen,GS_thresh,connectivity,nFirst=nMobile,iMax,relaxMethod)
%RelaxNetwork(): this function relaxes a node-filament network defined by
%the node coordinates (nodeXY), the filament connections (filNodeLUT), and the
%original filament lengths (oldFilLen) using relaxMethod relaxation method
%(currently only Gauss Seidel-esque method supported). Input parameter
%"connectivity" is determined before the start of stretching and assists in
%relaxation calculations.

%12/01/18: initial version split off from cellstretch_john_execute

%ERROR CHECKS


%/ERROR CHECKS


fSquare = 1;                        %represents summed squared forces (strain*E)^2 of all filaments
fSquareAll = zeros(iMax,1);         %stored fSquare values for each iteration

ii=1;                    %iteration counter, must be < iMax
strainFactor = 1;        %strainFactor is the factor that strain is multiplied by to convert it to step size
while (fSquare >= 0.1)   %ii++ every loop. Escape condition if ii > iMax
    
    %11/08/02
    [strcFilLen,strcFilAng]=FindFilLengthsAndAngles(nodeXY,filNodeLUT);   %find new filament lengths and angles
    
    %"Strain" is a vector of magnitudes for all filaments
    Strain = double((strcFilLen-oldFilLen)./oldFilLen);     %strcFilLen contains the "new" filament lengths. oldFilLen is the lengths before any stretch this cycle.
    
    %sum strains on the first nMobile nodes (the mobile nodes)
    nodeStrainXY = SumStrains(Strain,strcFilAng,connectivity);
    
    %adjust node positions (nodeXY) based on Strains
    [nodeXY,fSquare]=AdjustMovableNodes(nMobile,nodeXY,strainFactor,nodeStrainXY,elastMod);
    
    mobileNodeXY_new = mobileNodeXY_new + strainFactor.*nodeStrainXY;   %adjust node positions in the horizontal direction

    
    fSquareAll(ii)=fSquare;
    if ii >= iMax
        break
    elseif ii > 1
        strainFactor=AdjustStrainFactor(strainFactor,fSquareAll(ii-1),fSquareAll(ii));
    end
    ii = ii+1;  %iteration counter
end

if ii<iMax
    nodeXY_relaxed=nodeXY;  %this should be the relaxed network
else
    nodeXY_relaxed=[];  %network was not able to be relaxed
    msg=['WARNING (RelaxNetwork): network was not relaxed. iMax ' int2str(iMax) ' reached'];
end
