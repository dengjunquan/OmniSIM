function BSposAll = FindBSs(filename,centerLonLat)
str = fileread(filename);
pointstrings = loadjson(str);
numBSs = length(pointstrings.features);
BSxys  = cell(numBSs,1);

centerLongInDegrees = centerLonLat(1);
centerLatInDegrees = centerLonLat(2);

index = 0;

for i = 1:numBSs
    if ~isempty(pointstrings.features{i}.geometry)
        if strcmp(pointstrings.features{i}.geometry.type,'Point') == 1 && ...
           ~isempty(pointstrings.features{i}.geometry.coordinates(1))
            index = index+1;
            BSxys{index} = pointstrings.features{1,i}.geometry.coordinates;
        end
    end
end

numBSs = index;
BSxys(index+1:end)  = [];

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

BSposAll = zeros(numBSs,2);
for i = 1:numBSs
     BSxys{i}(:,1)  = (BSxys{i}(:,1) - centerLongInDegrees).*longlen;
     BSxys{i}(:,2)  = (BSxys{i}(:,2) - centerLatInDegrees).*latlen;
     BSposAll(i,:)  = BSxys{i};
end

end


