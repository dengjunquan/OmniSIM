function a = arrayresponse(theta,phy,obj)

M = obj.M;
N = obj.N;
dx = obj.dx;
dz = obj.dz;

option = 1;

if strcmp(obj.array,'UPA')==1
    if option
        u = exp((0:(M-1)).*dx.*1i.*2.*pi.*sin(theta));
        v = exp((0:(N-1)).*dz.*1i.*2.*pi.*sin(phy));
        a = kron(v,u);
        a = a.';          %////////////////////////////////////////////////
        
    else
        a = zeros(N*M,1);
        for m= 0:M-1
            for n= 0:N-1
                a(m*(N)+n+1) = exp(1i*2*pi*(m*dx*sin(phy)*cos(theta)+n*dz*cos(phy)));
            end
        end
    end
    
elseif strcmp(obj.array,'UCylinder')==1
    error('Uniform Cylinder Array will be implemented in the future!');
else
    error('unknown antenna array!');
end

end

