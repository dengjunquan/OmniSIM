function D = DiffractionCoeff(n,d_dash,d,phi_dash,phi,eta_0,eta_n,f)
c  = 299792458;
lamda = c/f;
k = 2*pi/lamda; % wavenumber
L = d_dash*d/(d_dash + d);
A = -exp(-1i*pi/4)/(2*n*sqrt(2*pi*k));
option = 1;
if option == 1 % see A New Heuristic UTD Diffraction Coefficient for Nonperfectly Conducting Wedges
    gamma1 = (pi-(phi-phi_dash))/(2*n);
    gamma2 = (pi+(phi-phi_dash))/(2*n);
    gamma3 = (pi-(phi+phi_dash))/(2*n);
    gamma4 = (pi+(phi+phi_dash))/(2*n);
    
    D1 = A*cot(gamma1)*TransitionF(2*k*L*n^2*(sin(gamma1))^2);
    D2 = A*cot(gamma2)*TransitionF(2*k*L*n^2*(sin(gamma2))^2);
    D3 = A*cot(gamma3)*TransitionF(2*k*L*n^2*(sin(gamma3))^2);
    D4 = A*cot(gamma4)*TransitionF(2*k*L*n^2*(sin(gamma4))^2);
end

if  option == 2  % see Cell Coverage Analysis of 28 GHz Millimeter Wave in Urban Microcell Environment Using 3-D Ray Tracing
    beta1 = phi-phi_dash;
    beta2 = phi+phi_dash;
    
    Nplus1 = round((pi+beta1)/(2*pi*n));
    Nplus2 = round((pi+beta2)/(2*pi*n));
    Nmius1 = round((-pi+beta1)/(2*pi*n));
    Nmius2 = round((-pi+beta2)/(2*pi*n));
    
    r1 = (pi+(phi-phi_dash))/(2*n);
    r2 = (pi-(phi-phi_dash))/(2*n);
    r3 = (pi+(phi+phi_dash))/(2*n);
    r4 = (pi-(phi+phi_dash))/(2*n);
    
    a1 = 2*(cos((2*n*pi*Nplus1 - beta1)/2))^2;
    D1 = A*cot(r1)*TransitionF(k*L*a1);
    a2 = 2*(cos((2*n*pi*Nmius1 - beta1)/2))^2;
    D2 = A*cot(r2)*TransitionF(k*L*a2);
    
    a3 = 2*(cos((2*n*pi*Nplus2 - beta2)/2))^2;
    D3 = A*cot(r3)*TransitionF(k*L*a3);
    a4 = 2*(cos((2*n*pi*Nmius2 - beta2)/2))^2;
    D4 = A*cot(r4)*TransitionF(k*L*a4);
end

% perpendicular polarization (parallel polarization)
theta_i = pi/2 - phi_dash; 
eta_r = eta_0;
R0 = (cos(theta_i) - sqrt(eta_r - sin(theta_i)^2))/(cos(theta_i)+ sqrt(eta_r - sin(theta_i)^2));
R0 = abs(R0);
theta_i = pi/2 - (n*pi-phi); 
eta_r = eta_n;
Rn = (cos(theta_i) - sqrt(eta_r - sin(theta_i)^2))/(cos(theta_i)+ sqrt(eta_r - sin(theta_i)^2));
Rn = abs(Rn);
D0 = D1+D2+R0*D3+Rn*D4;
D = abs(D0)*sqrt((d_dash+d)/(d_dash*d));
end

