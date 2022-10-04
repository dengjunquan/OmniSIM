function [shootpoint,wallid,phiout] = FindIntersection(startxy,phi,wallsMatrix)
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
raylength = 100000;
rayxys    = [startxy(1) + cos(phi)*1e-4       startxy(2) + sin(phi)*1e-4...
             startxy(1) + cos(phi)*raylength  startxy(2) + sin(phi)*raylength];

[Intersected, IntersectedXs, IntersectedYs] = lineIntersect(rayxys,wallsMatrix);

shootpoint = [NaN NaN];
wallid     = -1;
phiout     = NaN;

if sum(Intersected) > 0 % if there is a intersection point
    
    Distances = (startxy(1) - IntersectedXs).^2 + (startxy(2) - IntersectedYs).^2;
    Distances = Distances.*Intersected;
    
    Distances(Distances==0) = max(Distances);
    
    [~,wallids_sorted]  = sort(Distances);
    wallid = wallids_sorted(1);
    
    shootpoint = [IntersectedXs(wallid) IntersectedYs(wallid)];
    
    wallxys    = wallsMatrix(wallid,:);
    
    x1 = wallxys(1);
    y1 = wallxys(2);
    x2 = wallxys(3);
    y2 = wallxys(4);
    
    mirrorpoint = zeros(1,2);
    
    A = y2-y1;
    B = x1-x2;
    C = x1*(y1-y2) + y1*(x2-x1);
    
    a = startxy(1);
    b = startxy(2);
    
    mirrorpoint(1) = (a-2*A*(A*a+B*b+C)/(A*A+B*B));
    mirrorpoint(2) = (b-2*B*(A*a+B*b+C)/(A*A+B*B));
    
    [phiout,~] = cart2pol(shootpoint(1) - mirrorpoint(1), shootpoint(2) - mirrorpoint(2));
end
