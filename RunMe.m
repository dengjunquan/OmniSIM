close;clear;clc;
rng(888);
filename = 'Tsinghua.json';
BSxy = [-150 -100];

[walls,polygons] = FindWalls(filename);

wallsMatrix = zeros(length(walls),4);
wallsCenter = zeros(length(walls),2);
for walli = 1:length(walls)
    wallsMatrix(walli,1) = walls{walli,1}.p1(1);
    wallsMatrix(walli,2) = walls{walli,1}.p1(2);
    wallsMatrix(walli,3) = walls{walli,1}.p2(1);
    wallsMatrix(walli,4) = walls{walli,1}.p2(2);
    wallsCenter(walli,1) = (wallsMatrix(walli,1) + wallsMatrix(walli,3))/2;
    wallsCenter(walli,2) = (wallsMatrix(walli,2) + wallsMatrix(walli,3))/2;
end

ToPlot = 1;

if ToPlot ==1
    figure; hold on;
    % Draw walls
    if ~isempty(walls)
        for i=1:length(walls)
            plot([walls{i}.p1(1),walls{i}.p2(1)],[walls{i}.p1(2),walls{i}.p2(2)],'b');
        end
    end
    axis tight; axis equal;
    grid on;
end


cellRadius = 100;
numUEs = 100;
UEs = FindUEs(numUEs,polygons,cellRadius,BSxy);

testUEs = 1:numUEs;

if ToPlot == 1
    % Draw the UEs
    scatter(UEs(testUEs,1), UEs(testUEs,2), 2, jet(length(testUEs)),'filled');
end

% Draw BSs
if ToPlot ==1
    plot(BSxy(1),BSxy(2),'O'); hold on;
end

BouncingOrder = 8;
path_arrays = FastRT(walls,BSxy,UEs,BouncingOrder);


for UEindex = 8
   
    path_array = path_arrays{UEindex};
    for pathi = 1:length(path_array)
        path  = path_array{pathi};
        pathlen = size(path,1);
        for segmenti = 1:pathlen-1
            xy1 = path(segmenti,1:2);
            xy2 = path(segmenti+1,1:2);
            plot([xy1(1) xy2(1)],[xy1(2) xy2(2)],'-','Color',[0.1 0.1 0.1]);
        end
    end
end
plot(UEs(UEindex,1), UEs(UEindex,2), 'rO');
plot(BSxy(1),BSxy(2),'b^'); hold on;

