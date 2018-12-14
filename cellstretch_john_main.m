% This file contains functions:
%-cellstretch_john_main
%-ParametertCheck (does nothing right now)
%-GenerateReadme
%-UpdateReadme

function []=cellstretch_john_main(preFolder,datasetDescription,doElasticNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,stretchMag,breakMethod,breakFactor,intNodeNoise,nIntNodes,nPeriNodes,connectDensity,startRun,endRun,numCycles)
%12/01/06: changed the old parameter freeVerticalOffset to vertOffset since
%it no longer describes only which nodes are free.  Now, for the Stossel
%method, it describes the region that contains the nodes that are stretched
%and the nodes that are held stationary.  In the 11/12/18 update, this
%region defaulted to the upper 25% nodes (apical) and the bottom 25% nodes
%(basal) which would correspond to an offset of 45 deg. [PENDING].

%11/12/18: Two new parameter options.  Interior node positioning: 's' for square
%grid.  Stretch method: 's' for stretch method simulating Tom Stossell's
%Nature 2011 paper; this stretch method holds 25% of nodes (basal nodes)
%stationary and stretches 25% of nodes (apical nodes).

%11/09/06: Data is now saved in a subdirectory timestamped at the start of the
%main function

%11/09/05: takes 12 inputs now
%(datasetDescription,doElasticNodes,intNodePositioning,periNodePositioning,freeVer
%ticalOffset,breakMethod,nIntNodes,nPeriNodes,connectDensity,startRun,endRun,numCycles)

%11/08/03: takes 3 inputs now: nIntNodes, nPeriNodes, connectDensity

%11/07/27: need to check filament breakage algorithm as new filaments may
%not be replaced correctly

%-5/26/11: changed README formatting to coincide with new data organization
%techique.  Data will be put into folders using the e*_p*_o*_b* syntax
%where e describes node elasticity, p describes peripheral nodes, o
%describes vertical offset, and b describes break method employed.  Note
%that instead of denoting how many free vertical nodes there will be, it
%will be automatically calculated depending on the free vertical range.

%-1/20/11: added a doElasticNode option that does not reset positions of
%elastic nodes.  Also cleaned up functions including filament length and angle functions and the strain calculation function. Needs testing

%-1/13/11: noted that node positions are reset every cycle...will need to
%fix this to see if it makes a difference later.  Removed doSaveFig
%variable and instead use SAVE_DIR to make choice.  Also re-arranged when
%minicycles are recorded and saved to better match up with the figures.

%-12/28/10: finished cleaning up stretching and relaxing portions of code.  Data directories now created as "Data1", "Data2", etc. automatically to avoid overwriting

%-11/12/10: added initial support for "break low stress" method.  This breaks
%filaments experiencing low stress instead of vice versa.

%-11/8/10: bug...initial network is being overwritten every cycle

%-11/7/10: Added UpdateReadme function which updates data readme file when a
%run finishes and what cycle it ended on

%-11/5/10: added GenerateReadme function to output model parameters in same
%folder as data

%-11/1/10: broke up into a "main" and an "execute" function to decrease heap
%errors (hopefully)

%-9/22/10: will be used in resubmission to JTheorBiol
%-clear individual variables used after each run to prevent using incorrect
%data

%-8/30/10: will be used in resubmission to JTheorBiol
%-added a RandomizePeripheral function to create random but mirrored
%peripheral nodes
%-added a breakCounter variable to keep track of how many filaments have
%broken

%==Immediate To Do Changes==
%!!! Write high level pseudocode
%!!! Initial network is being overwritten every time!!! (FIXED)
%!Unify Fil->Node LUT to one LUT (NOT USEFUL because would have to search everytime)
%!Stop storing polar coordinates...store XY coordinates and convert to polar
%on the fly for graphing (DONE)
%!Track Angles (DONE)
%!Implement creation of internal nodes at intersections (DONE)

if nargin==0
    clear all; close all;
end

disp('Start cellstretch_john_main.m');
if ispc
    cd('C:\Users\John\Dropbox\PhD\Code\Cellstretch');
    data_dir='C:/Users/John/Dropbox/PhD/Data-Sim';
    dbstop if error %pauses execution if error is found
elseif isunix
    cd('~/cellstretch');
    data_dir='~/data-sim';
end

if nargin==0 %Default settings
    data_dir=pwd; %set current directory as data directory
    data_dir(data_dir=='\')='/'; %convert to Unix friendly format
    
    %Set model parameters
    preFolder = '141225 SR reviewer';  %subfolder that will contain further folders and subfolders.  Used to group experiments together
    
    datasetDescription=[preFolder ': altering # IN, PN, filaments'];   %will be inserted into readme file
    
    doElasticNodes=1; %if 1, ALL node positions are reset after every complete cycle.  If 0, then only fixed nodes are reset after every cycle (mobile nodes remain unmoved).
    
    intNodePositioning= 's'; %CHOICES: 'r' for random, 'i' for intersection of filaments, 's' for square grid
    
    periNodePositioning='u'; %CHOICES: 'r' for random, 'k' for kathy's positioning, 'u' for uniform
    
    %Describes the bidirectional offset from vertical (i.e., 270 and 90) on
    %the periphery.  For stretchMethod's other than Stossell's method, this
    %parameter will describe the region of free mobile nodes.  For
    %Stossel's method, this region holds the nodes that are stretched (top)
    %and the nodes that are held stationary (bottom).
    vertOffset=60;
    
    stretchMethod='s'; %CHOICES: 'p' for peripheral stretch/Kathy's method, 'i' for stretch subset of internal nodes and balance peripheral nodes (subsets determined by "vertOffset"), 's' for Stossel's shear method
    
    stretchMag=28; %Range is 0<integers<99: determines the magnitude of stretch for stretchMethod.  Default is 0.1.
    
    breakMethod='c';    %chooses method to break filaments. choices: kathy, threshold, none, break low stress, crit_conc
    
    breakFactor=1.0;    %factor that attenuates measures of breakage (individualized to different break methods)
    %NOTE: the higher the breakFactor (i.e. >1), the harder/less likely it is to break a filament
    
    %Set node and filament quantitaitive parameters
    intNodeNoise=0.03; % noise in the positioning of internal nodes. Represents fraction of cellRadius that is added as random noise
    
    %Note: set nIntNodes to [] to create maximum number of internal nodes
    %for a given nPeriNodes
    nIntNodes=[];    %number of nodes inside circle; these can move during stretch, =nnodes
    
    %number of peripheral nodes/focal adhesions. =nFocal
    nPeriNodes=40;  %should be divisible by 4 (for mirroring).  If periNodePositioning=='k', nPeriNodes must be 16.
    
    connectDensity=1.0;    %how well connected should the network be?  If 1, then maximally connected.
    
    %nPeriMobile = 6;  %this is UNKNOWN UNTIL later after peripheral nodes are setup; accounts for the nodes at the vertical regions that are mobile during peripheral stretch
    %filamentMultiplier = 15; %multiplies total # nodes by this to get # filaments
    
    %Initialize loop
    startRun=2;     %for labeling/data tracking purposes
    endRun = 5;   %for labeling/data tracking purposes
    numCycles = 1;   %how many times to stretch network stretchMag
    
end

%*****************************************
%******Initial model housekeeping******
%*****************************************

doRotate=true;  %rotates peripheral nodes by 15 degrees

%if nIntNodes not specified, will be set to the max number of internal
%nodes based on a uniform arrangement of peri nodes

%if nIntNodes not specified, set it to what it would be for a Stossel
%internal node network with uniform peripheral nodes (i.e., set to the max
%number of internal nodes based on a uniform arrangement of peri nodes)

if isempty(nIntNodes)
    if  stretchMethod=='s' && periNodePositioning=='u'
        nIntNodes=0;
        for ii=1:nPeriNodes/4 %for one hemisphere
            nIntNodes=nIntNodes+2*(ii-1);
        end
        nIntNodes=2*nIntNodes; %for other hemisphere
        
        if doRotate
            nIntNodes=nIntNodes+1;  %extra node is possible when peripheral nodes are rotated
        end
    elseif  stretchMethod=='s' && periNodePositioning=='s'
        nIntNodes = (nPeriNodes/4)^2;    
    end
end

nFilaments=CalcNumFilaments(intNodePositioning,nIntNodes,nPeriNodes,connectDensity,doRotate,periNodePositioning);

if nIntNodes==1 && nPeriNodes==4 && intNodePositioning=='s' %special case for "diamond" configuration
    nFilaments=8;
end

[doProceed, msg]=ParameterCheck(doElasticNodes, intNodePositioning, periNodePositioning, vertOffset, stretchMethod, stretchMag, breakMethod, breakFactor, doRotate, nIntNodes, nPeriNodes, connectDensity);

if ~doProceed
    disp('ERROR: ParameterCheck did not pass');
    disp(['Msg: ' msg]);
    return
end

%SETUP CONSTANTS
cellRadius = 100;   %radius of the undeformed Circle, =R
elastMod = 2e9;     %elastic modulus of filaments
delX = 1;           %amount of deformation applied to the circle each each cycle is stretchMag;
% relaxThresh=cellRadius/1e9;
relaxThresh=1e-3;
% totX = 10;          %amount of total deformation applied /REPLACED BY
% 'stretchMag' VARIABLE
iMax = 1e6;         %max number of interations for solutio n to converge

%doBiaxialStretch=0; %if 1, then vertical region nodes are mobile.  Otherwise, vertical regions nodes are fixed. LOOKUP WHETHER BV STRETCH IS UNIAXIAL OR BIAXIAL??

%===INITIALIZE SAVE/DISPLAY CONTROLS===
doDisplayFig=false; %whether figures are displayed...recommend to keep off for long runs
doSaveFig=true;   %whether figures are saved to file.  Recommend keep off for long runs
doSaveNetworkBeforeBreak=1; %determines whether *.mat file is saved before the filament break step

%   figOptions is two element vector.  First element controls whether
%   figure is displayed.  Second element controls whether figure is saved.
%       -to properly save figure, both figOptions(2) must be true and
%       SAVE_DIR must not be empty.  SAVEDIR=='./' saves in current
%       directory
figOptions=[doDisplayFig doSaveFig];

%===GENERATE FOLDER NAMES===
folderName=['e' int2str(doElasticNodes) '_i' intNodePositioning '_p' periNodePositioning '_o' int2str(vertOffset) '_s' stretchMethod(1) '_m' sprintf('%02u',stretchMag) '_b' breakMethod(1) '_f' sprintf('%05.2f',breakFactor) ];    %not used in naming
subFolderName=['IN' int2str(nIntNodes) '_PN' int2str(nPeriNodes) '_F' int2str(nFilaments) ];
start_time=GetCurrentTime();
ADDENDUM=['_intNodeNoise= ' num2str(intNodeNoise)];

if isempty(preFolder)
    BASE_DIR=[data_dir '/' folderName '/' subFolderName '/' start_time ADDENDUM];   %do not have to add '/' to end as CreateDir function will do that
else
    BASE_DIR=[data_dir '/' preFolder '/' folderName '/' subFolderName '/' start_time ADDENDUM];  %prefolder used to group experiments
end
[SAVE_DIR,exception]=CreateDir(BASE_DIR);

%Create ReadMe
GenerateReadme(SAVE_DIR,start_time,startRun,endRun,numCycles,datasetDescription,doElasticNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,stretchMag,breakMethod,breakFactor,doRotate,cellRadius,intNodeNoise,elastMod,nIntNodes,nPeriNodes,nFilaments,connectDensity,delX,relaxThresh,iMax,doSaveNetworkBeforeBreak);

tic_main=tic; %starts timer faor all runs
nTries=0;	%tracks how many attempts it took to execute a full run
whichRun=startRun;
while whichRun<=endRun %will keep looping until all runs finished
    tic;
    cellstretch_john_execute(SAVE_DIR,whichRun,numCycles,datasetDescription,doElasticNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,stretchMag,breakMethod,breakFactor,doRotate,cellRadius,intNodeNoise,elastMod,nIntNodes,nPeriNodes,nFilaments,delX,relaxThresh,iMax,doSaveNetworkBeforeBreak,figOptions);
    elapsedTime=toc;
    nTries=nTries+1;
    
    load([SAVE_DIR 'last_cycle']); %loads "lastCycle" variable
    
    if lastCycle>=numCycles  %if the execute script completed all cycles, update Readme and advance
        UpdateReadme(SAVE_DIR,start_time,elapsedTime,whichRun,lastCycle,nTries);
        whichRun=whichRun+1;
        nTries=0;
    end
end

%Note how long entire run lasted in ReadMe
elapsedTime_main=toc(tic_main);
elapsedHours_main=elapsedTime_main/3600;

existingReadme=fopen([SAVE_DIR 'data_readme_' start_time '.txt'],'a');      %append to end of file
fprintf(existingReadme,'Total runtime:  %3.1f hours \n',elapsedHours_main);
fclose(existingReadme);

clear all;
end

%ParameterCheck is currently a dummy function
function [doProceed, msg]=ParameterCheck(doElasticNodes, intNodePositioning, periNodePositioning, vertOffset, stretchMethod,stretchMag, breakMethod,breakFactor, doRotate, nIntNodes, nPeriNodes, connectDensity)

if ~exist('intNodePositioning','var') || isempty(intNodePositioning)
    
    switch intNodePositioning
        
        case 'i'
            
        case 's'
            nIntNodes_max=0;
            for ii=1:nPeriNodes/4 %for one hemisphere
                nIntNodes_max=nIntNodes_max+2*(ii-1);
                
                nIntNodes_max=2*nIntNodes_max; %for other hemisphere
            end
            if nIntNodes<nIntNodes_max
                msg='WARNING (cellstretch_main/ParameteCheck): Stossel network not maximally connected';
                disp(msg);
            end
    end
end

if ~exist('periNodePositioning','var') || isempty(periNodePositioning)
    
    switch periNodePositioning
        
    end
end

doProceed=true; %dummy function right now
msg=[];

end

function [] = GenerateReadme(SAVE_DIR,timestamp,startRun,endRun,numCycles,datasetDescription,doElasticNodes,intNodePositioning,periNodePositioning,vertOffset,stretchMethod,stretchMag,breakMethod,breakFactor,doRotate,cellRadius,intNodeNoise,elastMod,nIntNodes,nPeriNodes,nFilaments,connectDensity,delX,relaxThresh,iMax,doSaveNetworkBeforeBreak)
if nargin==0
    %SAVE_DIR=[];
    timestamp='1999-12-31_23;59';
    %timestamp='SAMPLE';
    datasetDescription='test readme';
    doElasticNodes=1;
    intNodePositioning='i';
    periNodePositioning='r';
    vertOffset=20;
    stretchMethod='stretch1';
    stretchMag=-10;
    breakMethod='break1';
    breakFactor=-121.521;
    startRun=1;
    endRun=10;
    numCycles=99;
    cellRadius=100;
    intNodeNoise=0.1;
    elastMod=1e9;
    nPeriNodes=10;
    nIntNodes=20;
    nFilaments=100;
    connectDensity=0.7;
    %     totX=10;  %REPLACED BY 'stretchMag' variable
    delX=1;
    relaxThresh=1e-99;
    iMax=2e9;
    doSaveNetworkBeforeBreak=false;
    %     SAVE_DIR= ['e' int2str(doElasticNodes) '_i' intNodePositioning '_p' periNodePositioning '_o' int2str(vertOffset) '_s' stretchMethod(1) '_m' sprintf('%02u',stretchMag) '_b' breakMethod(1) '_f' int2str(breakFactor) '/F' int2str(nFilaments) '_IN' int2str(nIntNodes) '_PN' int2str(nPeriNodes) '/' timestamp '/'];
    SAVE_DIR=['e' int2str(doElasticNodes) '_i' intNodePositioning '_p' periNodePositioning '_o' int2str(vertOffset) '_s' stretchMethod(1) '_m' sprintf('%02u',stretchMag) '_b' breakMethod(1) '_f' sprintf('%05.2f',breakFactor) '/F' int2str(nFilaments) '_IN' int2str(nIntNodes) '_PN' int2str(nPeriNodes) '/' timestamp '/' ];
    
    mkdir(SAVE_DIR);
end

%note whether sim run in Linux or Windows
if ispc
    OS='PC';
elseif isunix
    OS='Unix';
else
    OS='N/A';
end

%write readme file
newReadme=fopen([SAVE_DIR 'data_readme_' timestamp '.txt'],'w');
fprintf(newReadme,...
    [timestamp ...
    '\n' 'folder name:  %s'  ...
    '\n' 'dataset descripton: ' datasetDescription ...
    '\n' 'elastic nodes?: %0.1u'...
    '\n' 'internal node position type?: ', intNodePositioning,...
    '\n' 'peripheral node position type?: ',periNodePositioning,...
    '\n' 'vertical node offset: %02.0u' ...
    '\n' 'stretchMethod: ' stretchMethod ...
    '\n' 'stretchMag (%%): %02.0u' ...
    '\n' 'breakMethod: ' breakMethod ...
    '\n' 'breakFactor: %05.2f'...
    '\n' 'start run: %3.0u, putative end run: %3.0u' ...
    '\n' 'number of cycles: %3.0u'...
    '\n' 'cellRadius: %3.0u'...
    '\n' 'intNodeNoise: %3.3f' ...
    '\n' 'elastMod: %10.0u'...
    '\n' 'nIntNodes: %3.0u'...
    '\n' 'nPeriNodes: %3.0u'...
    '\n' 'nPeriMobile: depends on vertOffset and datasetType; see *.mat file'...
    '\n' 'nFilaments: %5.0u'...
    '\n' 'connectDensity: %0.4g'...
    '\n' 'delX: %2.0u'...
    '\n' 'relaxThresh: %g' ...
    '\n' 'iMax: %9.0u'...
    '\n' 'save network before filament breakage?: %1.0u'...
    '\n' 'OS: ' OS ...
    '\n']...
    ,SAVE_DIR,doElasticNodes,vertOffset,stretchMag,breakFactor,startRun,endRun,numCycles,cellRadius,intNodeNoise,elastMod,nIntNodes,nPeriNodes,nFilaments,connectDensity,delX,relaxThresh,iMax,doSaveNetworkBeforeBreak);

fclose(newReadme);
end

function []=UpdateReadme(SAVE_DIR,start_time,elapsedTime,lastRun,lastCycle,nTries)
%UpdateReadme: Reads the "last_cycle.mat" file saved in SAVE_DIR directory
%and inserts a line in the data readme indicating what cycle the last run
%(lastRun) terminated on

if nargin<5
    load([SAVE_DIR 'last_cycle']); %loads "last_cycle" variable
end

elapsedMinutes=elapsedTime/60;

existingReadme=fopen([SAVE_DIR 'data_readme_' start_time '.txt'],'a'); %append to end of file
timestamp=GetCurrentTime();
fprintf(existingReadme,['Run %2.0u finish on cycle %3.0u @ ' timestamp ' in %4.2f minutes after %d attempts \n'],lastRun,lastCycle,elapsedMinutes,nTries);
fclose(existingReadme);

end

%==Shortcomings==
%-only models focal adhesions on the periphery
%-doesn't take into account cooperative effect of microfilaments

%==Nodes/Focal adhesions==
%-"nPeriNodes/nFocal", "nIntNodes/nnodes" scalars represents # boundary focal adhesions (static) and nodes inside circle (mobile)
%-"periPol/afocal" 2xnPeriNodes holding polar cds of focal adhesions
%-"intPol/nodes" 2xnIntNodes holding polar cds of internal nodes
%-"nodeXY/XYall" is 2xnNodes matrix holding cartesian cds for nodes
%-"nodePol/nall" is 2xnNodes matrix holding polar cds for nodes

%==Filaments==
%-"filNodeLUT/filaments" 2xnFilaments, rep cxn between nodes; cols rep filament,
%row 1 values rep "origin" node, row 2 values rep "destination" node
%-"nodeFilLUT"/"nodefil" 36x36 rep filament connection (1-138) for each node; cols rep node; for each col/node, values rep filaments containing node
%-"numFils"/"lnodefil" 36x1 holds # nodes connected to each element
%-"filX"/"filY" 2x138 hold X,Y cd for filament end during stretch
%-"filR"/"filT" 2x138 hold polar cd for filament ends; these are generated at runtime typically from nodeXY and filNodeLUT at runtime
%-"filR0"/"filT0" are cds for initial filaments.
%-"originNode_numFils"/"scolumn" 36x1 with rows rep nodes and values rep # filaments with respective node as origin

%==Long term changes==
%-actin threshold for breakage in strain should be 20%(Janmey, JCB
%113:155-160,1991)
%-focal adhesions in center?-->no, also in periphery
%-actual radius of circle?-->need dimensions of fibroblasts
%-actual filament length? --> assume that 0.1-1uM long.  Contradicted by
%Protein Profile Actin
%strengthen filaments over time??

%==Questions
%*For code for breaking the strand, why use a random number generator to
%set threshold?
%*Where do the compressive forces come from in angles 80-100 and 260 to
%280?