function [walls,UEs,streets] = FindWallsandUEsNew(filename,numUEs,BS,cellRadius,UEspatialPat)

PositionRes = 5;

if strcmp(filename,'Manhattan')==1
    
    %First init the walls
    Nbx      = 6; Nby = 6;
    walls    = cell(Nbx*Nby*4,1);
    streets  = cell(6,1);
    polygons = cell(Nbx*Nby,1);
    
    index= 0;
    for xi = Nbx-1:-1:0
        for yi = Nby-1:-1:0
            
            Lby    = 50;  Lsx    = 40;
            Lbx    = 50;  Lsy    = 40;
            theta0 = unifrnd(-pi/24,pi/24);

            startX = 50+xi*(Lbx+Lsx)+unifrnd(-5,5);
            startY = 50+yi*(Lby+Lsy)+unifrnd(-5,5);
            
            W = Lby + unifrnd(-10,0); 
            L = Lbx + unifrnd(-10,0);
            
            p1 = [startX, startY];
            p2 = [startX, startY+W];
            p3 = [startX+L, startY+W];
            p4 = [startX+L, startY];
            
            if xi==Nbx-1 || xi ==0 
                theta0 = pi/2; % + unifrnd(-pi/24,pi/24);
            end
            
            ps = [p1;p2;p3;p4];
            center = sum(ps)/4;
            
            
            r = sqrt(L^2+W^2)/2;
            phy = atan(W/L);
            theta(1) = theta0 + phy;
            theta(2) = theta(1) + pi - 2*phy;
            theta(3) = theta(2) + 2*phy;
            theta(4) = theta(3) + pi - 2*phy;
            for j = 1:4
                ps(j,:) = [r*cos(theta(j))+ center(1),r*sin(theta(j))+ center(2)];
            end
            
            index=index+1;
            
            walls{index}.p1 = ps(1,:);
            walls{index}.p2 = ps(2,:);
            
            index=index+1;
            walls{index}.p1 = ps(2,:);
            walls{index}.p2 = ps(3,:);
            
            index=index+1;
            walls{index}.p1 = ps(3,:);
            walls{index}.p2 = ps(4,:);

            index=index+1;
            walls{index}.p1 = ps(4,:);
            walls{index}.p2 = ps(1,:);
            
            
            polygons{xi*Nby+yi+1}.xs = [ps(1,1);ps(2,1);ps(3,1);ps(4,1);ps(1,1)];
            polygons{xi*Nby+yi+1}.ys = [ps(1,2);ps(2,2);ps(3,2);ps(4,2);ps(1,2)];
        end
    end
    
    pp1 = [-(Lbx+Lsx) -(Lby+Lsy)] + 300;
    pp2 = [ (Lbx+Lsx) -(Lby+Lsy)] + 300;
    pp3 = [ (Lbx+Lsx)  (Lby+Lsy)] + 300;
    pp4 = [-(Lbx+Lsx)  (Lby+Lsy)] + 300;
    
    ppp1 = [-(Lbx+Lsx)+ 5             0] + 300;
    ppp2 = [0             -(Lbx+Lsx)+ 5] + 300;
    ppp3 = [(Lbx+Lsx)- 5              0] + 300;
    ppp4 = [0              (Lbx+Lsx)- 5] + 300;
    
    streets{1}.p1 = pp1; streets{1}.p1(1) = streets{1}.p1(1) + 5;
    streets{1}.p2 = pp2; streets{1}.p2(1) = streets{1}.p2(1) + 5;
    streets{2}.p1 = pp2; streets{2}.p1(2) = streets{2}.p1(2) + 5;
    streets{2}.p2 = pp3; streets{2}.p2(2) = streets{2}.p2(2) + 5;
    streets{3}.p1 = pp3; streets{3}.p1(1) = streets{3}.p1(1) - 5;
    streets{3}.p2 = pp4; streets{3}.p2(1) = streets{3}.p2(1) - 5;
    streets{4}.p1 = pp4; streets{4}.p1(2) = streets{4}.p1(2) - 5;
    streets{4}.p2 = pp1; streets{4}.p2(2) = streets{4}.p2(2) - 5;
    
    streets{5}.p1 = ppp1;
    streets{5}.p2 = ppp3;
    
    streets{6}.p1 = ppp2;
    streets{6}.p2 = ppp4;
    
    % UEs = zeros(numUEs,2);
    % Nbx = Nbx - 2;
    % Nby = Nby - 2;
    % for testi = 1:numUEs
    %     if rand()>0.5
    %         xl = rand()*(Lbx+Lsx)*Nbx;
    %         yl = mod(rand()*100-50 + randi(Nby,1)*(Lby+Lsy),Nby*(Lby+Lsy));
    %     else
    %         xl = mod(rand()*100-50 + randi(Nbx,1)*(Lbx+Lsx),Nbx*(Lbx+Lsx));
    %         yl = rand()*(Lby+Lsy)*Nby;
    %     end
    %     UEs(testi,:) = [xl+Lbx+Lsx,yl+Lby+Lsy];
    % end
        
elseif strcmp(filename,'RandomBuilding')==1
    density = 1.5e-04;
    Lmin = 40;
    Lmax = 80;
    Wmin = 20;
    Wmax = 40;
    Lx = 500;
    Ly = 500;
    Dmin = 30;
    [walls,polygons] = RandomBuilding(density,Lmin,Lmax,Wmin,Wmax,Lx,Ly,Dmin);
elseif strcmp(filename,'OpenSquare')==1
    walls = [];
    polygons = [];
    [UEs,~] = Uniformcircle(numUEs,10,cellRadius);
    UEs = UEs';
    UEs(:,1) = abs(UEs(:,1));
else
    
    index = 0;
    buildings = buildingFootprint(filename);
    walls = cell(length(buildings)*10,1);
    polygons = cell(length(buildings),1);
    for i = 1:length(buildings)
        polygon = buildings{i};
        xs = zeros(length(polygon),1);
        ys = zeros(length(polygon),1);
        for j = 1:length(polygon)-1
            index=index+1;
            walls{index}.p1 = polygon(j,:);
            walls{index}.p2 = polygon(j+1,:);
            xs(j) = polygon(j,1);
            ys(j) = polygon(j,2);
        end 
            xs(end) = polygon(1,1);
            ys(end) = polygon(1,2);
            polygons{i}.xs = xs;
            polygons{i}.ys = ys;
    end
    walls(index+1:end)=[];
end



if strcmp(UEspatialPat,'PPP')==1
    [UEs,~] = Uniformcircle(2*numUEs,10,cellRadius);
    UEs = UEs';
    UEs(:,1) = UEs(:,1) + BS(1);
    UEs(:,2) = UEs(:,2) + BS(2);
    
    UEsAll  = 1:2*numUEs;
    UEsOutdoor = UEsAll;
    for i = 1:length(polygons)
        [in,~]  = inpolygon(UEs(UEsOutdoor,1),UEs(UEsOutdoor,2),polygons{i}.xs,polygons{i}.ys);
        UEsOutdoor = setdiff(UEsOutdoor,UEsOutdoor(in));
    end
    UEsOutdoor = UEsOutdoor(1:numUEs);
    UEs        = UEs(UEsOutdoor,:);
    %UEs(:,1)   = abs(UEs(:,1));
    
elseif strcmp(UEspatialPat,'Grid') == 1
    
    Xs = -cellRadius:PositionRes:cellRadius;
    Ys = -cellRadius:PositionRes:cellRadius;
    XYs = allcomb(Xs,Ys);
    DistancesToBS = XYs(:,1).^2+XYs(:,2).^2;
    UEs = XYs(DistancesToBS<cellRadius.^2,:);
    
    UEs(:,1) = UEs(:,1) + BS(1);
    UEs(:,2) = UEs(:,2) + BS(2);
    
    UEsAll  = 1:length(UEs);
    UEsOutdoor = UEsAll;
    for i = 1:length(polygons)
        [in,~]  = inpolygon(UEs(UEsOutdoor,1),UEs(UEsOutdoor,2),polygons{i}.xs,polygons{i}.ys);
        UEsOutdoor = setdiff(UEsOutdoor,UEsOutdoor(in));
    end
    UEs        = UEs(UEsOutdoor,:);
    
elseif strcmp(UEspatialPat,'Trace') == 1
    
    if strcmp(filename,'Manhattan')==0
        error('1-D trace for Manhattan scenario only!');
    end
    Xs = [150 300 450];
    Ys = [125 225 325];
    
    UEs = [];
    
    Ystmp  = 125:PositionRes:325;
    Ystmp1 = 125:PositionRes:225;
    Ystmp2 = 226:PositionRes:325;
    
    Xstmp = 150:PositionRes:450;
    
    Y = 125;
    UEstmp = zeros(length(Xstmp),2);
    UEstmp(:,1) = Xstmp;
    UEstmp(:,2) = Y;
    UEs = [UEs;UEstmp];
    
    X = 450;
    UEstmp = zeros(length(Ystmp),2);
    UEstmp(:,1) = X;
    UEstmp(:,2) = Ystmp;
    UEs = [UEs;UEstmp];

    
%     X = 450;
%     UEstmp = zeros(length(Ystmp1),2);
%     UEstmp(:,1) = X;
%     UEstmp(:,2) = Ystmp1;
%     UEs = [UEs;UEstmp];
%     
%     Y = 225;
%     UEstmp = zeros(length(Xstmp),2);
%     UEstmp(:,1) = Xstmp(end:-1:1);
%     UEstmp(:,2) = Y;
%     UEs = [UEs;UEstmp];
%     
%     X = 150;
%     UEstmp = zeros(length(Ystmp2),2);
%     UEstmp(:,1) = X;
%     UEstmp(:,2) = Ystmp2;
%     UEs = [UEs;UEstmp];
%     
%     Y = 325;
%     UEstmp = zeros(length(Xstmp),2);
%     UEstmp(:,1) = Xstmp;
%     UEstmp(:,2) = Y;
%     UEs = [UEs;UEstmp];
    
elseif strcmp(UEspatialPat,'MaternCuster')==1
    AvgClusterSize = 5;
    r = 30;
    [Centers,~] = Uniformcircle(ceil(2*numUEs/AvgClusterSize),10,cellRadius);
    Centers = Centers';
    Centers(:,1) = Centers(:,1) + BS(1);
    Centers(:,2) = Centers(:,2) + BS(2);
    X = zeros(2*numUEs,2);
    total_so_far = 0;
        for c=1:length(Centers)
            NC = AvgClusterSize;           %poissrnd(AvgClusterSize); %number of points in cluster
            k = 0;
            while k < NC                   %draw uniformly in the n-ball via accept-reject
                Y = 2*r*rand(1,2) - r;     %candidate point
                if norm(Y) < r
                    X(total_so_far+k+1,:) = Centers(c,:) + Y;
                    k = k+1;
                end
            end
            total_so_far = total_so_far + NC;
        end
        UEs = X(1:total_so_far,:);         %cut off unused rows
        
        UEsAll  = 1:length(UEs);
        UEsOutdoor = UEsAll;
        for i = 1:length(polygons)
            [in,~]  = inpolygon(UEs(UEsOutdoor,1),UEs(UEsOutdoor,2),polygons{i}.xs,polygons{i}.ys);
            UEsOutdoor = setdiff(UEsOutdoor,UEsOutdoor(in));
        end
        UEsOutdoor = UEsOutdoor(1:numUEs);
        UEs        = UEs(UEsOutdoor,:);
            
elseif strcmp(UEspatialPat,'VANET')==1
    UEs=[];
    numUEsPerSteet = ceil(numUEs/6);
    streetwidth    = 10;
    
    moreUElength   = 0;
    
    for streetindex = 1:length(streets)
        p1 = streets{streetindex}.p1;
        p2 = streets{streetindex}.p2;
        
        if p1(1) == p2(1)
            streetlength   =  abs(p1(2) - p2(2)) + moreUElength*2;
            xys = [(rand(numUEsPerSteet,1)-1/2)*streetwidth + p1(1)  rand(numUEsPerSteet,1)*streetlength + min(p1(2),p2(2)) - moreUElength];
            %[~,SortedIndexs] = sort(xys(:,2)); 
            %xys = xys(SortedIndexs,:);
            UEs = [UEs; xys];
        elseif p1(2) == p2(2)
            streetlength   =  abs(p1(1) - p2(1)) + moreUElength*2;
            xys = [rand(numUEsPerSteet,1)*streetlength + min(p1(1),p2(1)) - moreUElength  (rand(numUEsPerSteet,1)-1/2)*streetwidth + p1(2)];
            %[~,SortedIndexs] = sort(xys(:,1)); 
            %xys = xys(SortedIndexs,:);
            UEs = [UEs; xys];
        end
    end
    UEs = UEs(1:numUEs,:);
        
end