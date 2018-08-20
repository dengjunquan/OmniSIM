function walls = buildingFootprint(jsonfilename)
% jsonfilename = 'manhattan.json';
str = fileread(jsonfilename);
buildings = loadjson(str);
numwalls = length(buildings.features);
walls = cell(numwalls,1);
index = 0;
for i = 1:numwalls
    if ~isempty(buildings.features{i}.geometry)
        if strcmp(buildings.features{i}.geometry.type,'Polygon') == 1 && ...
           ~isempty(buildings.features{i}.geometry.coordinates{1})
            index = index+1;
            walls{index} = buildings.features{i}.geometry.coordinates{1};
        end
    end
end

walls(index+1:end)  = [];

centerLongInDegrees = walls{round(length(walls))}(1,1);
centerLatInDegrees  = walls{round(length(walls))}(1,2);

% Set up "Constants"
m1 = 111132.92;     % latitude calculation term 1
m2 = -559.82;       % latitude calculation term 2
m3 = 1.175;         % latitude calculation term 3
m4 = -0.0023;       % latitude calculation term 4
p1 = 111412.84;     % longitude calculation term 1
p2 = -93.5;         % longitude calculation term 2
p3 = 0.118;         % longitude calculation term 3

% Calculate the length of a degree of latitude and longitude in meters
lat = centerLatInDegrees;
latlen = m1 + (m2 * cos(2 * lat)) + (m3 * cos(4 * lat)) + (m4 * cos(6 * lat));
longlen = (p1 * cos(lat)) + (p2 * cos(3 * lat)) + (p3 * cos(5 * lat));

for i = 1:length(walls)
     walls{i}(:,1)  = (walls{i}(:,1) - centerLongInDegrees).*longlen;
     walls{i}(:,2)  = (walls{i}(:,2) - centerLatInDegrees).*latlen;
end
end
