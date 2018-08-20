function UEs = FindUEs(numUEs,polygons,cellRadius,BS)
factor = 2;
[UEs,~] = Uniformcircle(factor*numUEs,10,cellRadius);
    UEs = UEs';
    UEs(:,1) = UEs(:,1) + BS(1);
    UEs(:,2) = UEs(:,2) + BS(2);
    
    UEsAll  = 1:2*numUEs;
    UEsOutdoor = UEsAll;
    for i = 1:length(polygons)
        [in,~]  = inpolygon(UEs(UEsOutdoor,1),UEs(UEsOutdoor,2),polygons{i}.xs,polygons{i}.ys);
        UEsOutdoor = setdiff(UEsOutdoor,UEsOutdoor(in));
    end
    if numUEs <= length(UEsOutdoor)
        UEsOutdoor = UEsOutdoor(1:numUEs);
        UEs        = UEs(UEsOutdoor,:);
    else
        error('try to increase the value of factor used in FindsUEs.m!');
    end
end

