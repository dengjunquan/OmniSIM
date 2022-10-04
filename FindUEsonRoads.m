function [UEs,TimeStamps] = FindUEsonRoads(Roads,UEdensity)
UEs = []; % Positions of UE Samples
TimeStamps = []; % TimeStamps of CSIs
%streetwidth    = 10;
t=-1000;
for streetindex = 1:length(Roads)
    p1 = Roads{streetindex}.p1;
    p2 = Roads{streetindex}.p2;
    
    L = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2);
    N = ceil(L/UEdensity);
    %R = rand(1,N);
    R = 0:(1/N):1;
    xys = [p1(1)+R.*(p2(1)-p1(1)); p1(2)+R.*(p2(2)-p1(2))];  
    
    UEs = [UEs xys];

    if streetindex>1 && p1(1)==Roads{streetindex-1}.p2(1) && p1(2)==Roads{streetindex-1}.p2(2)
        TimeStamps=[TimeStamps (t+1):t+length(xys)];
    else
        t = t+1000;
        TimeStamps=[TimeStamps (t+1):t+length(xys)];
    end
    t= t+length(xys);
end
end

