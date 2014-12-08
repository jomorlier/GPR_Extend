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
mean_ = {@meanZero};
hyp0.mean = [];
lik = {@likGauss};
sn = 0.2;
hyp0.lik = log(sn);
inf = {@infExact};





% hold on;
% %gradient descent produces strange results
% for i = 1:1000
% [nlZ, dnlZ] = gp(hyp0, inf, mean_, cov, lik, xtrain, ytrain);
% disp(dnlZ)
% hyp0.cov = hyp0.cov-0.0001*dnlZ.cov;
% hyp0.lik = hyp0.lik-0.0001*dnlZ.lik;
% disp(hyp0)
% 
% disp(mean(mean(nlZ)))
% plot(i,mean(mean(nlZ)))
% plotyy(i,dnlZ.cov(2),'g')
% end


