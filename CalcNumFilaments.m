%called by:
%-GenerateNetwork

function [nFilaments]=CalcNumFilaments(intNodePositioning,nIntNodes,nPeriNodes,connectDensity,doRotate,periNodePositioning)
%CalcNumFilaments: calculates the expected number of total filaments given
%the type of internal node positioning, the # of internal nodes, the # of
%peripheral nodes, and how densely connected the graph will be.

%12/01/16: add support for parameter doRotate to accurately calculate the
%correct number of filaments

%11/12/27: add support for Stossel method...shares same code as the
%intersection method

if nargin==0
    intNodePositioning='s';
        periNodePositioning='s';

    nIntNodes=30;
    nPeriNodes=80;
    if intNodePositioning=='s'
        nIntNodes=0;
        for ii=1:nPeriNodes/4 %for one hemisphere
            nIntNodes=nIntNodes+2*(ii-1);
        end
        nIntNodes=2*nIntNodes; %for other hemisphere
    end
    connectDensity=138/nchoosek(nPeriNodes+nIntNodes,2);
end

switch intNodePositioning
    case 'r' %random
        %nFilaments=filamentMultiplier*(nPeriNodes+nIntNodes); old method
        nFilaments=uint32(nchoosek(nPeriNodes+nIntNodes,2)*connectDensity); %filaments generated randomly interconnecting any two nodes (peri or internal)
        
    case {'i'} %intx b/w peripheral-node-to-peripheral-node filaments
        %Total # filaments takes into account the number of ways to connect
        %nPeriNodes (nchoosek term) AND how many filaments will be
        %generated via intersections (2*nIntNodes term)
        nFilaments=uint32(nchoosek(nPeriNodes,2)*connectDensity+2*nIntNodes);  %See lab notebook for details.
        
    case {'s'} %square grid (Stossel 2011 method)
        %Calculates filaments in one hemisphere in one direction (i.e. top
        %hemisphere, horizontal direction) and multiples by four to get
        %other hemispheres.  Adds the filaments b/w adjacent nodes at the
        %end
        nFilaments=0;
        
        if periNodePositioning == 's' %flat top and bottom per SR reviewer
            
            for ii=1:nPeriNodes/4
                nFilaments=nFilaments+ii;
            end
            nFilaments=nFilaments*4;    %reflect to other three quadrants
            nFilaments=nFilaments+nPeriNodes;   %connect filaments on periphery
            nFilaments=nFilaments*connectDensity;
            
        else
            if doRotate     %rotated peripheral node arrangement by 15 degrees
                for ii=1:nPeriNodes/4-1
                    nFilaments=nFilaments+2*ii;
                end
                nFilaments=nFilaments*4;    %reflect to other three quadrants
                nFilaments=nFilaments+nPeriNodes;   %connect filaments on periphery
                nFilaments=nFilaments+2*2*nPeriNodes/4; %filaments on the equators
                nFilaments=nFilaments*connectDensity;
            else    %original peripheral node arrangement
                
                for ii=2:nPeriNodes/4
                    nFilaments=nFilaments+2*ii-1;   %for the grid filaments in first hemisphere (e.g., top)
                end
                nFilaments=nFilaments*4;        %reflect to other three hemispheres (e.g., bot, left, right)
                nFilaments=nFilaments+nPeriNodes; %for the filaments b/w adjacent peri nodes
                nFilaments=nFilaments*connectDensity;
            end
            
        end
        
        
        
        
end

end