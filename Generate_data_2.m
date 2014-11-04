%Generate data 2
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

cov = {@covSumBasis, {@covSEiso, @covConst_2}}; sf = 1; ell = 0.4;                             % setup the GP
% hyp0.cov  = [log([ell;sf]) log([ell;sf])]; 
c2 = 4;
sf2 = 2;
B = 1;
hyp0.cov = [log([ell;sf]) ;log(c2); sf2; B];
% hyp0.cov = [];
mean = {@meanSum, {@meanZero,@meanZero}};
hyp0.mean = [];
lik = {@likGauss}; sn = 0.2;
hyp0.lik  = log(sn);
inf = {@infExact}; 
Ncg = 50;                                   % number of conjugate gradient steps
% ymu{1} = f(xte); ys2{1} = sn^2; nlZ(1) = -Inf;

hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtrain', ytrain'); % opt hypers
[ymu ys2 fmu fs2] = gp(hyp, inf, mean, cov, lik, xtrain', ytrain', xtest');

hold on;
plot(xtrain,ytrain,'o','MarkerSize',8);
% plot(xtest,ytest,'o');
plot(xtest,ymu,'.r');
plot(xtest,ytest,'.c')
[ymu ys2 fm2 fs2 lp] = gp(hyp, inf, mean, cov, lik, xtrain', ytrain', xtest', ytest');
xlabel('Explanatory');
ylabel('Output');
title('Comparison of Actual Function, Model Output');
legend('Train','Model Output','Underlying Function');