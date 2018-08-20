function d = DistanceToLineSegment(P0,P1,P2)
d01 = sqrt( (P0(1)-P1(1))^2 +  (P0(2)-P1(2))^2 );
d02 = sqrt( (P0(1)-P2(1))^2 +  (P0(2)-P2(2))^2 );
d12 = sqrt( (P2(1)-P1(1))^2 +  (P2(2)-P1(2))^2 );

x1 = P1(1);
y1 = P1(2);
x2 = P2(1);
y2 = P2(2);

A = y2-y1;
B = x1-x2;
C = x1*(y1-y2) + y1*(x2-x1);

d = abs(A*P0(1) + B*P0(2) + C)./sqrt(A^2+B^2);

if d01^2 >= d02^2+d12^2
    d = d02;
end

if d02^2 >= d01^2+d12^2
    d = d01;
end
end

