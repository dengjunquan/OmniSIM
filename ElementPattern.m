function G = ElementPattern(theta,phy,obj)
if strcmp(obj.ElementPattern,'omnidirectional') ==1
    G = 1;
elseif strcmp(obj.ElementPattern,'patchcosine') ==1
    G = sqrt(5.3428)*cos(theta)*cos(phy);
    if  cos(theta) < 0
        G = 0;
    end
elseif strcmp(obj.ElementPattern,'sector') ==1
    G = 2;
    if  cos(theta) < 1/4
        G = 0;
    end
else
    error('unknow antenna element pattern!');
end

% =============== Test =================
% MinG = 1/3;
% thetas = -pi:pi/180:pi;
% phys   = -pi:pi/180:pi;
% [Thetas,Phys] = meshgrid(thetas,phys);
% xs = sin(Phys).*cos(Thetas);
% ys = sin(Phys).*sin(Thetas);
% zs = cos(Phys);
% 
% Alphas = acos(ys./(sqrt(xs.^2+ys.^2+zs.^2)));
% ElePats = max(MinG,abs(cos(Alphas))).^2;
% figure;
% mesh(Thetas,Phys,ElePats);
% figure;
% mesh(Thetas, Phys, max(MinG,abs(sin(Thetas).*sin(Phys)).^2));

end

