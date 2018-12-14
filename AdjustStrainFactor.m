function [strainFactor_new] = AdjustStrainFactor(strainFactor,prevVal,currVal)
%AdjustStrainFactor
%   Adjusts the strainFactor (controls the amount of change in node positions)

%   If strainFactor is not small and stresses are increasing, strainFactor
%   is decreased.

%   If strainFactor is not large and stresses are increasing OR decrease
%   too slowly, strainFactor is increased.

%12/01/23: changed input parameters fSquare_prev and fSquare_current to
%prevVal and currVal respectively.

if (strainFactor >0.1) && (currVal > prevVal)	%if strain factor not small AND stresses inc, then decrease strain factor
    strainFactor_new = 9/10*strainFactor;
    %RUSSELL: SUGGEST TO TRY TO REVERT TO PREVIOUS STATE --> DON'T THINK
    %THIS WOULD WORK BECAUSE THERE IS NO STOCHASTICITY AFTER NETWORK
    %GENERATION UP TO THIS POINT. STOCASTICITY IS IN THE BREAKAGE.
    
elseif (strainFactor <2.5) && (1.005*currVal > prevVal)	%if strain factor not large AND stresses inc or dec too sowly, increase strain factor
    strainFactor_new = 10/9*strainFactor;
else
    strainFactor_new=strainFactor;
end
end