function [LUT_new] = GenerateRandomFilaments(nFil,nAll,LUT_old)
%GenerateRnadomFilaments:
%   -Generate specified number (nFil) of random, unique filaments (LUT_new)
%   given # of filaments (nFil) and # of nodes (nAll) and optionally, an
%   old table of filaments (LUT_old).
%
%   -if LUT_old is given as input, nFil is the number of new
%   filaments desired and nAll is the total number of nodes available, and
%   LUT_new is the look up table for the new filaments only
%
%   -NOTE: there are checks for filament duplication (i.e. two filaments
%   may have the same node pairs) and for self-filaments (i.e. same node
%   connected to itself)

%   Code used for 1) initial generation of filaments in GenerateNetwork, 2)
%   generation of replacement filaments after breakage in
%   cellstretch_execute

%11/12/27: name change from GenerateFilaments to GenerateRandomFilaments

%11/08/01: added optional 3rd parameter LUT_old such that new filaments
%can be generated that are not already in LUT_old

%INPUT CHECKS
if nargin==0
    clear all; close all;
    disp('Start Generate Filaments');
    nFil=3;
    nAll=5;
    LUT_old=[1 2; 2 3; 1 3]';
elseif nargin<3
    LUT_old=[];
end

%ERROR CHECKS
if nFil==0 || nAll==0
    disp('ERROR (GenerateFilaments): number of filaments desired or number of nodes is 0');
    LUT_new=[];
    return
elseif nFil>nchoosek(nAll,2)-size(LUT_old,2) %max # filaments = all combinations of any two nodes
    disp(['ERROR (GenerateFilaments): desired ',int2str(nFil),' filaments but only ',int2str(nchoosek(nAll,2)-size(LUT_old,2)),' possible']);
    LUT_new=[];
    return
end

nRemain=nFil;    %# new fils to be created. Is altered later.
LUT_new=[];
while nRemain>0
    
    %Create nRemain new filaments
    newFils=NaN(2,nRemain);
    newFils(1,:) = unidrnd(nAll,1,nRemain);   %random sequence of all nodes for "origin" nodes
    newFils(2,:) = unidrnd(nAll,1,nRemain);   %random sequence of all nodes for "destination" node
    
    %Concatenate new filament
    LUT_new=[LUT_new,newFils];

    %****Check that there are no self-filaments (i.e. same node connected to itself)****
    for ii = 1:size(LUT_new,2)
        while LUT_new(1,ii)==LUT_new(2,ii)
            LUT_new(2,ii) = unidrnd(nAll,1,1);
        end
    end
    
    %****Check that filaments aren't duplicated****
    LUT_new=sort(LUT_new);  %-ascending sort every column so that union can properly be calculated
    LUT_old=sort(LUT_old);  %/
    
    if ~isempty(LUT_old) 
        LUT_uni=union(LUT_new',LUT_old','rows');
    else
        LUT_uni=union(LUT_new',LUT_new','rows'); %union of cols of matrix in LUT_new (i.e. removes duplicate cols). Sorted.
    end
    LUT_new=LUT_uni';

    nRemain=nFil+size(LUT_old,2)-size(LUT_new,2); %# filaments that remain to be created.
    %Note that size of LUT_u is >= size of LUT_old
 
end

end
