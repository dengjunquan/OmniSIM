function DiffEdgesFinal = FindDiffractionEdges(BSxy,wallsMatrix)
%raylength = 100000;

ToPlot = 0;
numWalls = length(wallsMatrix);

phis = (0:0.005:360-0.00001).*pi./180;  
numPhis = length(phis);
IntersectedPoints  = zeros(numPhis,2);
IntersectedWalls  =  zeros(numPhis,1);

for phiIndex = 1:numPhis
    phi = phis(phiIndex);
    startxy = BSxy;
    [shootpoint,wallid,~] = FindIntersection(startxy,phi,wallsMatrix);
     if wallid > 0
            IntersectedPoints(phiIndex,:)  = shootpoint;
            IntersectedWalls(phiIndex)     = wallid;
     end
end

wallids = unique(IntersectedWalls); % LOS Walls

if wallids(1) == 0
    wallids = wallids(2:end);
end

% 判断 IntersectedPoints 前后是否连续
Pe = IntersectedPoints(end,:);
P1 = IntersectedPoints(1,:);
P2 = IntersectedPoints(2,:);

IndexMark = NaN;
if  norm(Pe-P1) < 1.8 * norm(P1-P2) % 前后连接
    for ii = 1:numPhis-2
        P1 = IntersectedPoints(ii,:);
        P2 = IntersectedPoints(ii+1,:);
        P3 = IntersectedPoints(ii+2,:);
        if norm(P3-P2) > 1.8 * norm(P1-P2)
            IndexMark = ii+1;
        break;
        end
    end
end

%disp(IndexMark);
IntersectedPoints = cat(1,IntersectedPoints(IndexMark+1:end,:),IntersectedPoints(1:IndexMark,:));
IntersectedWalls  = cat(1,IntersectedWalls(IndexMark+1:end,:), IntersectedWalls(1:IndexMark,:));

%figure;hold on;
LoSwallsMatrix = zeros(length(wallids),4);
for i = 1:length(wallids)
    wallid = wallids(i);
    %SegPointsonWall = IntersectedPoints(IntersectedWalls == wallid,:);
    %disp(SegPointsonWall);
    LoSwallsMatrix(i,:) = wallsMatrix(wallid,:);
    %disp(wallsMatrix(wallid,:));
    %plot(wallsMatrix(wallid,1:2:3),wallsMatrix(wallid,2:2:4),'linewidth',2);
    %text(wallsMatrix(wallid,1),wallsMatrix(wallid,2),num2str(wallid));
end

if ToPlot == 1
    for i=1:length(wallsMatrix)
        plot(wallsMatrix(i,1:2:3),wallsMatrix(i,2:2:4),'linewidth',1);
    end
    
    hold on;
    scatter(IntersectedPoints(:,1),IntersectedPoints(:,2),'.');
    hold on;
    scatter(BSxy(1),BSxy(2),'O');
end

DiffEdges = [];
for ii = 1:numPhis-3
    P1 = IntersectedPoints(ii,:);
    P2 = IntersectedPoints(ii+1,:);
    P3 = IntersectedPoints(ii+2,:);
    P4 = IntersectedPoints(ii+3,:);
    
    d1 = norm(P2-P1);
    d2 = norm(P3-P2);
    d3 = norm(P4-P3);
    ds = sort([d1 d2 d3]);
    
    if norm(P3-P2) > 1.8 * ds(2)  % 找到不连续点 
        Edge = [P2 IntersectedWalls(ii+1) d1];
        DiffEdges = cat(1,DiffEdges,Edge);
        Edge = [P3 IntersectedWalls(ii+2) d3];
        DiffEdges = cat(1,DiffEdges,Edge);
    end
end

d1 = norm(IntersectedPoints(2,:) - IntersectedPoints(1,:));
Edge = [IntersectedPoints(1,:) IntersectedWalls(1) d1];
DiffEdges = cat(1,DiffEdges,Edge);

d3 = norm(IntersectedPoints(end,:) - IntersectedPoints(end-1,:));
Edge = [IntersectedPoints(end,:) IntersectedWalls(end) d3];
DiffEdges = cat(1,DiffEdges,Edge);


DiffEdgesFinal = DiffEdges;
index = 1;
for jj = 1:length(DiffEdges(:,1))
    Edge = DiffEdges(jj,:);
    if Edge(3)==0
        continue;
    end
    wall = wallsMatrix(Edge(3),:);
    p = Edge(1:2);
    d = Edge(4);
    wpoint1 = wall(1:2);
    wpoint2 = wall(3:4);
    d1 = norm(p-wpoint1);
    d2 = norm(p-wpoint2);
    dmin = min([d1 d2]);
    if dmin == d1
        diffpoint = wpoint1;
    else
        diffpoint = wpoint2;
    end
    % find n face
    zeroface = Edge(3);
    nface = 0;
    for wi = 1:numWalls
       wall = wallsMatrix(wi,:);
       wpoint1 = wall(1:2);
       wpoint2 = wall(3:4);
       if (norm(wpoint1-diffpoint)==0||norm(wpoint2-diffpoint)==0) && wi~= Edge(3)
           nface = wi;
       end
    end
    if dmin < 1.8*d % Find a diffraction edge
        diffedge = [diffpoint zeroface nface];
        DiffEdgesFinal(index,:) = diffedge;
        index = index + 1;
    end
end
DiffEdgesFinal(index:end,:) = [];
if ToPlot == 1
    scatter(DiffEdgesFinal(:,1),DiffEdgesFinal(:,2),'*');
end
