function [xy,rt]=Uniformcircle(n,r1,r2)
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
% Set your seed here if desired
%  n = 50000;
%  r1 = 5; r2 = 10;
 r = sqrt(r1^2+(r2^2-r1^2)*rand(1,n)); % Using square root here ensures distribution uniformity by area
 t = 2*pi*rand(1,n);
 x = r.*cos(t);
 y = r.*sin(t);
 xy=[x;y];
 rt=[r1*r,t];
%  plot(x,y,'.');
end