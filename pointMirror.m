function mirrorpoint =  pointMirror(p0,p1,p2)

x1 = p1(1);
y1 = p1(2);
x2 = p2(1);
y2 = p2(2);

mirrorpoint = zeros(1,2);

A = y2-y1;
B = x1-x2;
C = x1*(y1-y2) + y1*(x2-x1);

a = p0(1);
b = p0(2);

mirrorpoint(1) = (a-2*A*(A*a+B*b+C)/(A*A+B*B));
mirrorpoint(2) = (b-2*B*(A*a+B*b+C)/(A*A+B*B));

% p0 = [0 0]; p1 = [1 2]; p2 = [3 2];
% pm = pointMirror(p0,p1,p2);
% figure; hold on;
% plot(p0(1),p0(2),'o'); plot(pm(1),pm(2),'o');
% plot([ p1(1) p2(1)], [p1(2) p2(2)],'-'); axis equal;

end