%This file contains functions:
%1. FindIntersections_Schwarz
%2. intersections

function [intxXY,intxLines,nIntx_final]=FindIntersections_Schwarz(nIntx_desired,ptsXY,Lines)
%FindIntersections: given a number of intersections to be found (nIntx_desired), the
%coordinates of the points (ptsXY) and how points are connected
%(Lines), returns the XY-crds of the intersections (intxXY) and
%the IDs of the lines involved in the intersections (intxLines).
%
%NOTE: intersections AT the points themselves are not counted, i.e., the
%intersection must not be at an endpoint of a line
%
%11/07/11: this version uses "intersections" by Paul Schwarz.  It gives
%correct intersections.  Note that nIntx_desired is not guaranteed and thus
%function returns nIntx_final<nIntx_desired if there are insufficient #
%intx

if nargin==0
    close all; clear all;
    example=3;
    switch example
        case(1)
            nIntx_desired=2;    %number of intersections to look for
            ptsXY=[0.7 -0.7 0.7 -0.7 1 0; 0.7 0.7 -0.7 -0.7 0 1];
            Lines=[1 2 5 6; 4 3 4 1];    %opposite vertices of square are connected
        case(2)
            nIntx_desired=2;    %number of intersections to look for
            ptsXY=[70.7106781186548,-70.7106781186547,-70.7106781186548,70.7106781186547;70.7106781186547,70.7106781186548,-70.7106781186547,-70.7106781186548];
            Lines=[1,1,1,2,2,3;2,3,4,3,4,4];
        case (3)
            nIntx_desired=1;
            ptsXY=[70.7106781186548,-70.7106781186547,-70.7106781186548,70.7106781186547;70.7106781186547,70.7106781186548,-70.7106781186547,-70.7106781186548];
            Lines=[1 2; 3 4]';
    end
    plot(ptsXY(1,:),ptsXY(2,:),'*','color','red');    %plots points as red stars
    hold on;
    for ii=1:size(Lines,2)  %plot lines
        plot(ptsXY(1,Lines(:,ii)),ptsXY(2,Lines(:,ii)));
        text(linspace(ptsXY(1,Lines(1,ii)),ptsXY(1,Lines(2,ii)),10),linspace(ptsXY(2,Lines(1,ii)),ptsXY(2,Lines(2,ii)),10),int2str(ii),'Color','black','FontWeight','normal'); %label lines with line IDs
        hold on;
    end
    xLim=[1.5*min(min(ptsXY)) 1.5*max(max(ptsXY))]; yLim=[1.5*min(min(ptsXY)) 1.5*max(max(ptsXY))];
    set(gca,'XLim',xLim,'YLim',yLim);
end

%INPUT CHECKS
intxXY=[]; intxLines=[]; nIntx_final=[];
if max(max(Lines))>size(ptsXY,2) || min(min(Lines))<1
    disp('ERROR (FindIntersections): lines refer to unknown points');
    return
elseif nIntx_desired<1
    disp(['WARNING (FindIntersections): n=',int2str(nIntx_desired),' intersection points desired']);
    nIntx_final=0;
    return
elseif size(Lines,2)<2
    disp(['WARNING (FindIntersections): n=',int2str(size(Lines,2)),' lines cannot have any intersections']);
    return
elseif size(ptsXY,2)<3
    disp(['WARNING (FindIntersections): n=',int2str(size(ptsXY,2)),' points cannot have any intersections']);
    return
elseif nIntx_desired>nchoosek(size(Lines,2),2)
    disp(['WARNING (FindIntersections): ',int2str(nIntx_desired),' intersections cannot be generated from ',int2str(size(Lines,2)),' lines']);
end

% n=1e3; %Part of troubleshooting loop
% foundIntx=zeros(n,1); %Part of troubleshooting loop
% numTrials=zeros(n,1); %Part of troubleshooting loop
nIntx_final=nIntx_desired; %may change below if not enough intx found

% for ii=1:n %Part of troubleshooting loop
intxXY=NaN(2,nIntx_desired);    %stores coordinates of intersection. row1=X-crds, row2=Y-crds
intxLines=NaN(2,nIntx_desired);    %stores lines that are intersecting. row1=line1, row2=point2's
nLines=size(Lines,2); %num of lines

trial=1;
iMax=1e3*factorial(size(ptsXY,2))/factorial(size(ptsXY,2)-2);    %# permutations of 2 points from all possible pts x10
intx_idx=1;

%Is this truly a random intersection finder?
while intx_idx<=nIntx_desired
    line1=0; line2=0;
    while (line1==line2) %choose two non-equivalent lines at random
        line1=randi(nLines); line2=randi(nLines);
    end
    l1_XY=[ptsXY(:,Lines(1,line1)) ptsXY(:,Lines(2,line1))]; %-each line is defined by 4 coordinates: 2 X- and 2 Y-crds
    l2_XY=[ptsXY(:,Lines(1,line2)) ptsXY(:,Lines(2,line2))]; %/
    
    if nargin==0; disp(['Trial #',int2str(trial),'. l1=',int2str(line1),'. l2=',int2str(line2)]); end;
    
    [intxX,intxY]=intersections(l1_XY(1,:),l1_XY(2,:),l2_XY(1,:),l2_XY(2,:),1); %if no intx found, [intxX, intxY]=empty
    
    isAnIntx=unique(~isempty([intxX;intxY]));  %isAnIntx is true if an intersection was found by func "intersections".
    
    %isNewIntxPt is positive if both conditions true:
    %   (1) intx crds are not already stored in intxXY
    %   (2) intx crds are not already stored in ptsXY (i.e., intx @ existing pt). 
    %Note: is true when [intxX,intxY] is empty
    isNewIntxPt=isempty(find(FindColumnInMatrix([ptsXY intxXY],[intxX;intxY]),1));
    
    if isAnIntx && isNewIntxPt     %if lineintersect finds an intersection and this intersection is new....
        if nargin==0; disp(['Found! ',int2str([intxX intxY])]); end
        intxXY(:,intx_idx)=[intxX;intxY];    %....store the intersection crds and...
        intxLines(:,intx_idx)=sort([line1;line2]);         %...stores the lines IDs that have intersected.  Ascending sorted.
        intx_idx=intx_idx+1;  %...and increment the counter
        %numTrials(ii)=trial; %Part of troubleshooting loop
    end
    trial=trial+1;
    if trial>iMax
        nIntx_final=sum(~isnan(mean(intxXY)));  %gives final number of intersections stored in intxXY
        disp(['WARNING (FindIntersections): not able to find ',int2str(nIntx_desired),' intersections. Only ',int2str(nIntx_final), ' intersections found.  Try increasing iMax']);
        break;
    end
end

%     foundIntx(ii)=sum(~isnan(mean(intxXY))); %Part of troubleshooting loop
% end %Part of troubleshooting loop

%Remove NaNs from the ends (on right) of intxXY and intxLines
[~,c]=find(isnan(intxXY));
intxXY(:,c)=[];
intxLines(:,c)=[];

if nargin==0
    plot(intxXY(1,:),intxXY(2,:),'o','MarkerFaceColor','g','LineWidth',2);  %plots intersection with blue 'o'
end
% sum(foundIntx==1) %Part of troubleshooting loop
% avgNumTrials=mean(numTrials) %Part of troubleshooting loop
end

function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.

% Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.


% Input checks.
error(nargchk(2,5,nargin))

% Adjustments when fewer than five arguments are supplied.
switch nargin
	case 2
		robust = true;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 3
		robust = x2;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 4
		robust = true;
		self_intersect = false;
	case 5
		self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
		length(x1) ~= length(y1)
	error('X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
		length(x2) ~= length(y2)
	error('X2 and Y2 must be equal-length vectors of at least 2 points.')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
	repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
	repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
	repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
	repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
	repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Force i and j to be column vectors, even when their length is zero, i.e.,
% we want them to be 0-by-1 instead of 0-by-0.
i = reshape(i,[],1);
j = reshape(j,[],1);

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.

% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.

if robust
	overlap = false(n,1);
	warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
	try
		lastwarn('')
		for k = 1:n
			T(:,k) = AA(:,:,k)\B(:,k);
			[unused,last_warn] = lastwarn;
			lastwarn('')
			if strcmp(last_warn,'MATLAB:singularMatrix')
				% Force in_range(k) to be false.
				T(1,k) = NaN;
				% Determine if these segments overlap or are just parallel.
				overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			end
		end
		warning(warning_state)
	catch err
		warning(warning_state)
		rethrow(err)
	end
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
	% For overlapping segment pairs the algorithm will return an
	% intersection point that is at the center of the overlapping region.
	if any(overlap)
		ia = i(overlap);
		ja = j(overlap);
		% set x0 and y0 to middle of overlapping region.
		T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		selected = in_range | overlap;
	else
		selected = in_range;
	end
	xy0 = T(3:4,selected).';
	
	% Remove duplicate intersection points.
	[xy0,index] = unique(xy0,'rows');
	x0 = xy0(:,1);
	y0 = xy0(:,2);
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		sel_index = find(selected);
		sel = sel_index(index);
		iout = i(sel) + T(1,sel).';
		jout = j(sel) + T(2,sel).';
	end
else % non-robust option
	for k = 1:n
		[L,U] = lu(AA(:,:,k));
		T(:,k) = U\(L\B(:,k));
	end
	
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		iout = i(in_range) + T(1,in_range).';
		jout = j(in_range) + T(2,in_range).';
	end
end

% Plot the results (useful for debugging).
% plot(x1,y1,x2,y2,x0,y0,'ok');
end