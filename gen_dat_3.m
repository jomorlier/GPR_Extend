%Gen_dat_3
clear;
x = linspace(0,30,1000);
xpost = linspace(30,50,1000);
r = randperm(length(x));
x = x(r);
ntrain = 20;
y = cos([x xpost])+-0.5*sqrt([x xpost]);
xtrain = x(1:ntrain);
xtest = [x(ntrain+1:end) xpost];
ytrain = y(1:ntrain);
ytest = y(ntrain+1:end);

cov = {@covSEiso};
sf =  1;
ell = 0.4;
hyp0.cov = [log(sf), log(ell)];
mean = {@meanZero};
hyp0.mean = [];
lik = {@likGauss};
sn = 0.2;
hyp0.lik = log(sn);
inf = {@infExact};

% disp(hyp0.cov)
% [nlZ, dnlZ] = gp(hyp0, inf, mean, cov, lik, xtrain, ytrain);
% disp(dnlZ)
% hyp0.cov = hyp0.cov-0.1*dnlZ.cov;
% disp(hyp0.cov)
% disp(sum(sum(nlZ)))
% [nlZ, dnlZ] = gp(hyp0, inf, mean, cov, lik, xtrain, ytrain);
% disp(dnlZ)
% hyp0.cov = hyp0.cov-0.00001*dnlZ.cov;
% disp(hyp0.cov)
% disp(sum(sum(nlZ)))
% [nlZ, dnlZ] = gp(hyp0, inf, mean, cov, lik, xtrain, ytrain);
% disp(dnlZ)
% disp(sum(sum(nlZ)))


