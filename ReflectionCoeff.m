function coeff = ReflectionCoeff(theta_i,eta_r)
%==========================================================================
% 1- Schlick's approximation
% eta_r = 5;
% n1 = 1.0003;
% n2 = sqrt(eta_r);
% R0 =((n1-n2)/(n1+n2 ))^2;
% Rs = zeros(1,90);
% for theta = 1:90
%     Rs(theta) = R0+(1-R0 )*(1-cos(theta.*pi/180)).^5;
% end
% plot(1:90,Rs);

%  2 - https://en.wikipedia.org/wiki/Fresnel_equations
%  1) 73 GHz Wideband Millimeter-Wave Foliage and Ground Reflection Measurements and Models .
%  2) A Comparison of Theoretical and Empirical Reflection Coefficients for Typical Exterior  Wall Surfaces in a Mobile Radio Environment.
%  3) Time-Domain Free-Field Measurements of the Relative Permittivity of Building Materials.
%  4) Measurements of Reflection and Transmission Characteristics of Interior Structures of Office Building in the 60GHz Band.
%  n1*sin(theta_i) = n2*sin(theta_t)  (1)
%  Cs = (n1*cos(theta_i)-n2*cos(theta_t))/(n1*cos(theta_i)+n2*cos(theta_t))
%     = (cos(theta_i)-sqrt(eta_r - sin(theta_i)^2))/(cos(theta_i)-sqrt(eta_r + sin(theta_i)^2))
%  Cp = (n1*cos(theta_t)-n2*cos(theta_i))/(n1*cos(theta_t)+n2*cos(theta_i))
%     = (sqrt(eta_r-sin(theta_i)^2) - eta_r*cos(theta_i))/(sqrt(eta_r-sin(theta_i)^2) + eta_r*cos(theta_i))
%  pho_s = exp(-8*(pi*delta_h*cos(theta_i)/lamda)^2);
%  Csavg  = ((1+pho_s)/2)*Cs;
%  Cpavg  = ((1+pho_s)/2)*Cp;
%  Cavg   = ((1+pho_s)/2)*(Cs+Cp)/2;

% delta_h = 0.0025;
% c  = 299792458;
% fc = 28e9;
% lamda = c/fc;
% coeffs=zeros(2,5,90);
% eta_rs = 4:8;
% for j=1:5
%     eta_r = eta_rs(j);
%     for i = 1:90
%         theta_i = i*pi/180;
%         Cs = -(cos(theta_i) - sqrt(eta_r - sin(theta_i)^2))/(cos(theta_i)+ sqrt(eta_r - sin(theta_i)^2));
%         Cp =  (sqrt(eta_r-sin(theta_i)^2) - eta_r*cos(theta_i))/(sqrt(eta_r-sin(theta_i)^2) + eta_r*cos(theta_i));
%         pho_s  = exp(-8*(pi*delta_h*cos(theta_i)/lamda)^2);
%         % coeff  = ((1+pho_s)/2)*(abs(Cs)+abs(Cp))/2;
%         coeff  = ((1+pho_s)/2)*abs(Cs);
%         coeffs(1,j,i)=coeff;
%         coeff  = ((1+pho_s)/2)*abs(Cp);
%         coeffs(2,j,i)=coeff;
%     end
% end
% figure;hold on;
% for j=1:5
%     plot(1:90,squeeze(coeffs(1,j,:)));
%     plot(1:90,squeeze(coeffs(2,j,:)));
% end

%==========================================================================
delta_h = 0.0025;
c  = 299792458;
fc = 28e9; 
lamda = c/fc;
Cs = -(cos(theta_i) - sqrt(eta_r - sin(theta_i)^2))/(cos(theta_i)+ sqrt(eta_r - sin(theta_i)^2));
% Cp =  (sqrt(eta_r-sin(theta_i)^2) - eta_r*cos(theta_i))/(sqrt(eta_r-sin(theta_i)^2) + eta_r*cos(theta_i));
pho_s  = exp(-8*(pi*delta_h*cos(theta_i)/lamda)^2);
% coeff  = ((1+pho_s)/2)*(abs(Cs)+abs(Cp))/2;
% coeff  = ((1+pho_s)/2)*abs(Cp);
coeff  = ((1+pho_s)/2)*abs(Cs);
end