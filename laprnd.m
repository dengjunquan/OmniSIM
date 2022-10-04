function y = laprnd(n, mu, sigma, max)
% generation of a numbers with truncated Laplace distribution
index=0;
y = zeros(1,n);
if sigma ==0
    y = mu+y;
    return;
end
while index<n
    x = log(rand()./rand());
    if abs(sigma*x)<max
        index=index+1;
        y(index) = mu + sigma*x;
    end
end