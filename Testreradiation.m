alpha = 0.8;
beta = pi*20/180;
thetas = -pi:0.001:pi;
P =zeros(1,length(thetas));
for i = 1:length(thetas)
    theta = thetas(i);
    P(i)  = alpha*(2/beta)^2*exp(-(theta/beta)^2) + 1-alpha;
end
P = P./(alpha*(2/beta)^2 + 1-alpha);
figure;
plot(180.*thetas./pi,10.*log10(P));