function f = TransitionF(X)
% A Uniform Geometrical Theory of Diffraction for an Edge in a Perfectly Conducting Surface
% % when X is small
% f1 = (sqrt(pi*X)-2*X*exp(1i*pi/4)-(2/3)*X^2*exp(-1i*pi/4))*exp(1i*(pi/4 + X));
% % when X is large
% f2 = 1+1i*(1/(2*X))-(3/4)*(1/X^2)-1i*(15/8)*(1/X^3)+(75/16)*(1/X^4) - 1i*(375/32)*(1/X^4);

% Xs = [0.001:0.001:0.01 0.02:0.01:0.1 0.2:0.1:1 2:1:10 20:10:100];
% Fs1 = zeros(length(Xs),1);
% Fs2 = zeros(length(Xs),1);
% Fs =  zeros(length(Xs),1);
% for ii = 1:length(Xs)
%     X = Xs(ii);
%     Fs1(ii) = (sqrt(pi*X)-2*X*exp(1i*pi/4)-(2/3)*X^2*exp(-1i*pi/4))*exp(1i*(pi/4 + X));
%     Fs2(ii) = 1+1i*(1/(2*X))-(3/4)*(1/X^2)-1i*(15/8)*(1/X^3)+(75/16)*(1/X^4) - 1i*(375/32)*(1/X^4);
%     Fs(ii) = 2*1i*sqrt(X)*exp(1i*X)*integral(@(t)exp(-1i*t.^2),sqrt(X),400);
% end
% figure;
% semilogx(Xs,abs(Fs1),'-');hold on;grid on;
% semilogx(Xs,abs(Fs2),'-');
% semilogx(Xs,abs(Fs),'-');
% ylim([0 1]);
% figure;
% semilogx(Xs,angle(Fs1),'-');hold on;grid on;
% semilogx(Xs,angle(Fs2),'-');
% semilogx(Xs,angle(Fs),'-');
% ylim([0 1]);

if X<=0.3
    f = (sqrt(pi*X)-2*X*exp(1i*pi/4)-(2/3)*X^2*exp(-1i*pi/4))*exp(1i*(pi/4 + X));
elseif X>=6
    f = 1+1i*(1/(2*X))-(3/4)*(1/X^2)-1i*(15/8)*(1/X^3)+(75/16)*(1/X^4) - 1i*(375/32)*(1/X^4);
else
    f = 2*1i*sqrt(X)*exp(1i*X)*integral(@(t)exp(-1i*t.^2),sqrt(X),280);
end

