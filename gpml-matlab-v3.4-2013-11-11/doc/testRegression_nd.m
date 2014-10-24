addpath('../../..'); %to add basis functions in higher folder

disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')

disp(' '); disp('clear all, close all')
clear all, close all
write_fig = 0;
disp(' ')

%Define the mean, covariance, and likelihood functions
disp('meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];')
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hypml = [1;1;1]; hyp.mean = [hypml; 1];
disp('covfunc = {@covMaternard, 3}; ellp = 1/4; ellv = 1/4; ellu = 1/4; sf = 1; hyp.cov = log([ellp; ellv; ellu; sf]);')
covfunc = {@covMaternard, 3}; ellp = 1/4; ellv = 1/4; ellu = 1/4; sf = 1; hyp.cov = log([ellp; ellv; ellu; sf]);
disp('likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);')
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn); 
disp(' ')

%Generate dataset input and underlying function parameters
[sim_ds,E] = orig_sim;%orig_sim(0.01,100,pi/12,1,1.05,0,pi/2,0);
xprev = sim_ds(1:(end/30),2:end);
timeset = sim_ds(1:(end/30),1);
r = randperm(length(xprev(:,1)));
x = xprev(r,:);
timesetr = timeset(r);
n = length(xprev(:,1));
l = n/12;                                %size of training set
%freq = 12;                             %frequency for underlying function
%amp = .75;                             %amplitude for underlying function

traininds = [1:l];                      %indicies for different parts of data
testinds = [l+1:n];
testinds = [traininds testinds];        %to test over the whole distribution

%Define the underlying function and output w/ noise
%disp('x = gpml_randn(0.3, n, 1);')
%x = gpml_randn(0.3, n, 1);
disp('K = feval(covfunc{:}, hyp.cov, x);')
K = feval(covfunc{:}, hyp.cov, x(traininds,:));
disp('mu = feval(meanfunc{:}, hyp.mean, x);')
mu = feval(meanfunc{:}, hyp.mean, x(traininds,:));
disp('ytr = chol(K)''*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);')
ftr = diff(xprev(:,1))';                         %outputs are derivatives
ftr = ftr(1,r(traininds));
noise = randn(1,length(ftr));                      %gaussian noise
ytr = ftr' + .005*noise';                         %training output w/ noise

%Plot dataset and underlying function
figure(2)
set(gca, 'FontSize', 24)
disp(' '); disp('plot(x(traininds), ytr, ''+'')')
hold on;
plot(timesetr(traininds,1), ytr(traininds,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
axis([min(timeset(:,1)) max(timeset(:,1)) -.1 .1])
grid on
xlabel('time, seconds')
ylabel('output, ytr')
if write_fig, print -depsc f1.eps; end
disp(' ');

%Generate gaussian process w/ training data
disp(' ')
disp('nlml = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr)')
nlml = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds,:), ytr)
disp(' ')

%Generate prediction w/ training data and test input
disp('z = xprev;')
z = xprev;
disp('[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds,:), ytr, z);

%Plot initial prediction
figure(3)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];') 
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8);')
fill([timeset(:,1); flipdim(timeset(:,1),1)], f, [7 7 7]/8);

disp('hold on; plot(z, m); plot(x, ytr, ''+'')')
hold on; plot(timeset(:,1), m, 'LineWidth', 2);
plot(timesetr(traininds,1), ytr(:,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
axis([min(timeset(:,1)) max(timeset(:,1)) -.1 .1])
grid on
xlabel('time, seconds')
ylabel('output, ytr')
if write_fig, print -depsc f2.eps; end
disp(' ');

disp(' ')
disp('covfunc = @covSEiso; hyp2.cov = [0; 0]; hyp2.lik = log(0.1);')
covfunc = @covSEard; hyp2.cov = [1; 1; 1; 1]; hyp2.lik = log(0.1);
disp('hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x, ytr)')
hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x(traininds,:), ytr);
disp(' ')

disp('exp(hyp2.lik)')
exp(hyp2.lik)
disp('nlml2 = gp(hyp2, @infExact, [], covfunc, likfunc, x, ytr)')
nlml2 = gp(hyp2, @infExact, [], covfunc, likfunc, x(traininds,:), ytr)
disp('[m s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x(traininds,:), ytr, z);

%Optimizing the hyperparameters
disp(' ')
figure(4)
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)');
fill([timeset(:,1); flipdim(timeset(:,1),1)], f, [7 7 7]/8)
disp('hold on; plot(z, m); plot(x, ytr, ''+'')');
hold on; plot(timeset(:,1), m, 'LineWidth', 2);

plot(timesetr(traininds,1), ytr(:,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('time, seconds')
ylabel('output, ytr')
axis([min(timeset(:,1)) max(timeset(:,1)) -.1 .1])
if write_fig, print -depsc f3.eps; end
disp(' ');

disp(' ')
disp('hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);')
hyp.cov = [1; 1; 1; 1]; hyp.mean = [[1;1;1]; 1]; hyp.lik = log(0.1);
disp('hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, ytr);')
hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x(traininds,:), ytr);
disp('[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds,:), ytr, z);

% Refine the prediction
figure(5)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];')
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)')
fill([timeset(:,1); flipdim(timeset(:,1),1)], f, [7 7 7]/8)
disp('hold on; plot(z, m); plot(x, ytr, ''+'');')
hold on; plot(timeset(:,1), m, 'LineWidth', 2);

plot(timesetr(traininds,1), ytr(:,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('time, seconds')
ylabel('output, ytr')
axis([min(timeset(:,1)) max(timeset(:,1)) -.1 .1])
if write_fig, print -depsc f4.eps; end
disp(' ');

disp('large scale regression using the FITC approximation')
disp('nu = fix(n/2); u = linspace(.6,3.4,nu)'';')
nu = fix(n/2); pu = linspace(.6,3.4,nu)'; vu = linspace(.6,3.4,nu)'; uu = linspace(.6,3.4,nu)'; u = [pu,vu,uu]; 
disp('covfuncF = {@covFITC, {covfunc}, u};')
covfuncF = {@covFITC, {covfunc}, u(traininds,:)};
disp('[mF s2F] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x, ytr, z);')
[mF s2F] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x(traininds,:), ytr, z);

%Adjust training bounds
figure(6)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [mF+2*sqrt(s2F); flipdim(mF-2*sqrt(s2F),1)];')
f = [mF+2*sqrt(s2F); flipdim(mF-2*sqrt(s2F),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)')
fill([timeset(:,1); flipdim(timeset(:,1),1)], f, [7 7 7]/8)
disp('hold on; plot(z, mF); plot(x, ytr, ''+'');')
hold on; plot(timeset(:,1), mF, 'LineWidth', 2);

plot(timesetr(traininds,1), ytr(:,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
disp('plot(u,1,''o'')')
plot(pu,1,'ko', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('time, seconds')
ylabel('output, ytr')
axis([min(timeset(:,1)) max(timeset(:,1)) -.1 .1])
if write_fig, print -depsc f5.eps; end
disp(' ')