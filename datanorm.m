function [y, Hs] = datanorm(x)
h = size(x,1);
mu = mean(x, 2);
sigma2 = sum(var(x, 0, 2));
s = sqrt((h-1)/sigma2);
Hs = [s*[eye(h-1), -s*mu(1:h-1)]; [zeros(1,h-1) 1]];
y = Hs*x;