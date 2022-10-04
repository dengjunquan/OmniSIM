clc;
% Xs = [0.001:0.001:0.01 0.02:0.01:0.1 0.2:0.1:5 5.2:0.2:100];
% Fs =  zeros(length(Xs),1);
% for ii = 1:length(Xs)
%     X = Xs(ii);
%     Fs(ii) = TransitionF(X);
% end
% figure;
% semilogx(Xs,abs(Fs),'-','linewidth',2);hold on;grid on;
% ylim([0 1]);
% 
% figure;
% semilogx(Xs,angle(Fs),'-','linewidth',2);hold on;grid on;
% ylim([0 1]);

etas = 3:6;
Neta = length(etas);
Phis = pi:0.01:(pi+pi/2);
N = length(Phis);
Ds = zeros(Neta,N);

n = 3/2;
phi_dash = pi/4;
f = 2.4e9;
d_dash = 10;
d = 10;
figure;hold on;grid on;
for ei = 1:Neta
    eta_0 = etas(ei);
    eta_n = etas(ei);
    for i = 1:N
        phi = Phis(i);
        Ds(ei,i) = DiffractionCoeff(n,d_dash,d,phi_dash,phi,eta_0,eta_n,f);
    end
    plot(Phis,Ds(ei,:));
end