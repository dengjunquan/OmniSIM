function path_arrays = FastRT(walls,BSxy,UEs,BouncingOrder)
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
phis = (0:0.1:360-0.00001).*pi./180;
numPhis = length(phis);
numUEs  = length(UEs);

wallsMatrix = zeros(length(walls),4);
wallsCenter = zeros(length(walls),2);
for walli = 1:length(walls)
    wallsMatrix(walli,1) = walls{walli,1}.p1(1);
    wallsMatrix(walli,2) = walls{walli,1}.p1(2);
    wallsMatrix(walli,3) = walls{walli,1}.p2(1);
    wallsMatrix(walli,4) = walls{walli,1}.p2(2);
    wallsCenter(walli,1) = (wallsMatrix(walli,1) + wallsMatrix(walli,3))/2;
    wallsCenter(walli,2) = (wallsMatrix(walli,2) + wallsMatrix(walli,4))/2;
end

IntersectedPoints  = zeros(BouncingOrder,numPhis,2);
BouncingDirections = zeros(BouncingOrder,numPhis);
IntersectedWalls  =  zeros(BouncingOrder,numPhis);

for phiIndex = 1:numPhis
    for  bouncingIndex = 1:BouncingOrder
        if bouncingIndex == 1
            phi = phis(phiIndex);
            startxy = BSxy;
        else
            phi = BouncingDirections(bouncingIndex-1,phiIndex);
            startxy = IntersectedPoints(bouncingIndex-1,phiIndex,:);
        end
        
        [shootpoint,wallid,phiout] = FindIntersection(startxy,phi,wallsMatrix);
        
        if wallid > 0
            IntersectedPoints(bouncingIndex,phiIndex,:)  = shootpoint;
            BouncingDirections(bouncingIndex,phiIndex)   = phiout;
            IntersectedWalls(bouncingIndex,phiIndex)     = wallid;
        end
        
        if wallid == -1
            raylength  = 1000;
            shootpoint = [startxy(1) + cos(phi)*raylength  startxy(2) + sin(phi)*raylength];
            IntersectedPoints(bouncingIndex,phiIndex,:)  = shootpoint;
            BouncingDirections(bouncingIndex,phiIndex)   = phiout;
            IntersectedWalls(bouncingIndex,phiIndex)     = wallid;
            break;
            
        end
        
    end
end

ToPlot = 0;
if ToPlot==1
for phiIndex = 1:100:numPhis
    for bouncingIndex = 1:BouncingOrder
        if bouncingIndex == 1
            intpointxy = IntersectedPoints(bouncingIndex,phiIndex,:);
            plot([BSxy(1) intpointxy(1)], [BSxy(2) intpointxy(2)],'-');
        elseif IntersectedWalls(bouncingIndex,phiIndex) > 0 
            intpointxy_pre = IntersectedPoints(bouncingIndex-1,phiIndex,:);
            intpointxy = IntersectedPoints(bouncingIndex,phiIndex,:);
            plot([intpointxy_pre(1) intpointxy(1)], [intpointxy_pre(2) intpointxy(2)],'-');
        elseif IntersectedWalls(bouncingIndex,phiIndex) == -1
            intpointxy_pre = IntersectedPoints(bouncingIndex-1,phiIndex,:);
            intpointxy = IntersectedPoints(bouncingIndex,phiIndex,:);
            plot([intpointxy_pre(1) intpointxy(1)], [intpointxy_pre(2) intpointxy(2)],'-');
        end
    end
end
end

path_arrays = cell(numUEs,1);

for UEindex = 1:numUEs
    UExy = UEs(UEindex,:);
    numpaths = 0;
    path_array = cell(100,1);
    path_IDs   = zeros(100,1);
    pathi      = 0;
    pathID     = 0;
    
    wallsOnPath_pre = []; % IDs for the walls which are on a ray path
    
    for phiIndex = 1:numPhis
        dxys = ones(BouncingOrder,1).*1000;
        for bouncingIndex = 1:BouncingOrder
            if bouncingIndex == 1
                intpointxy = IntersectedPoints(bouncingIndex,phiIndex,:);
                ray = [BSxy(1)  BSxy(2) intpointxy(1) intpointxy(2)];
            elseif IntersectedWalls(bouncingIndex,phiIndex) > 0
                intpointxy_pre = IntersectedPoints(bouncingIndex-1,phiIndex,:);
                intpointxy = IntersectedPoints(bouncingIndex,phiIndex,:);
                ray = [intpointxy_pre(1) intpointxy_pre(2) intpointxy(1) intpointxy(2)];
            elseif IntersectedWalls(bouncingIndex,phiIndex) == -1
                intpointxy_pre = IntersectedPoints(bouncingIndex-1,phiIndex,:);
                intpointxy = IntersectedPoints(bouncingIndex,phiIndex,:);
                ray = [intpointxy_pre(1) intpointxy_pre(2) intpointxy(1) intpointxy(2)];
            end
           
            dxys(bouncingIndex) = DistanceToLineSegment(UExy, ray(1:2), ray(3:4));
            
        end
        [dmin,bounce_order] = min(dxys);
        
        wallsOnPath_now = zeros(bounce_order + 1,1);
        
        if dmin < 1
            pathi = pathi + 1;
            numpaths = numpaths + 1; % Found a path
            path = zeros(bounce_order + 1,3);
            path(1,:) = [BSxy';0];
            wallsOnPath_now(1) = 0;
            for iorder = 1:bounce_order
                if iorder == bounce_order
                    path(iorder+1,:) = [UExy'; IntersectedWalls(iorder,phiIndex)];
                    wallsOnPath_now(iorder+1) = IntersectedWalls(iorder,phiIndex);
                else
                    path(iorder+1,:) = [squeeze(IntersectedPoints(iorder,phiIndex,:)); IntersectedWalls(iorder,phiIndex)];
                    wallsOnPath_now(iorder+1) = IntersectedWalls(iorder,phiIndex);
                end
            end
            if ~isequal(wallsOnPath_now,wallsOnPath_pre)
                pathID = pathID + 1; % a new path found!
            end
            wallsOnPath_pre   = wallsOnPath_now;
            
            path_array{pathi} = path;
            path_IDs(pathi)   = pathID;
        end
    end
    path_array(pathi+1:end) = [];
    path_IDs(pathi+1:end)   = [];
    
    path_array_merged = cell(pathID,1);
    
    for ii = 1:pathID
        indexs = find(path_IDs == ii);
        for jj=1:length(indexs)
            if jj==1
                path_array_merged{ii} = path_array{indexs(jj)};
            else
                path_array_merged{ii} = path_array_merged{ii} +  path_array{indexs(jj)};
            end
        end
        path_array_merged{ii} = path_array_merged{ii}./length(indexs);
    end
    
    path_arrays{UEindex} = path_array_merged;
end
end