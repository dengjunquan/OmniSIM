function [walls,polygons] = FindWalls(filename)
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
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

