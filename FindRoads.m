function  roadsegment= FindRoads(filename,centerLonLat)
str = fileread(filename);
linestrings = loadjson(str);
numroads = length(linestrings.features);
roads = cell(numroads,1);

centerLongInDegrees = centerLonLat(1);
centerLatInDegrees = centerLonLat(2);

index = 0;

for i = 1:numroads
    if ~isempty(linestrings.features{i}.geometry)
        if strcmp(linestrings.features{i}.geometry.type,'LineString') == 1 && ...
           ~isempty(linestrings.features{i}.geometry.coordinates(1))
            index = index+1;
            roads{index} = linestrings.features{1,i}.geometry.coordinates;
        end
    end
end

roads(index+1:end)  = [];

% Set up "Constants"
m1 = 111132.92;     % latitude calculation term 1
m2 = -559.82;       % latitude calculation term 2
m3 = 1.175;         % latitude calculation term 3
m4 = -0.0023;       % latitude calculation term 4
p1 = 111412.84;     % longitude calculation term 1
p2 = -93.5;         % longitude calculation term 2
p3 = 0.118;         % longitude calculation term 3

% Calculate the length of a degree of latitude and longitude in meters
lat = centerLatInDegrees.*pi/180;
latlen = m1 + (m2 * cos(2 * lat)) + (m3 * cos(4 * lat)) + (m4 * cos(6 * lat));
longlen = (p1 * cos(lat)) + (p2 * cos(3 * lat)) + (p3 * cos(5 * lat));

for i = 1:length(roads)
     roads{i}(:,1)  = (roads{i}(:,1) - centerLongInDegrees).*longlen;
     roads{i}(:,2)  = (roads{i}(:,2) - centerLatInDegrees).*latlen;
end

roadsegment = cell(length(roads)*20,1);

index = 0;
for i = 1:length(roads)
    road = roads{i};
    for j = 1:length(road)-1
        index=index+1;
        roadsegment{index}.p1 = road(j,:);
        roadsegment{index}.p2 = road(j+1,:);
    end
end
roadsegment(index+1:end)=[];
end

