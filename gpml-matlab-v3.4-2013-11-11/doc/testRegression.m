disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
disp('Hit any key to continue...'); pause

disp(' '); disp('clear all, close all')
clear all, close all
write_fig = 0;
disp(' ')

%Define the mean, covariance, and likelihood functions
disp('meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];')
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
disp('covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);')
covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
disp('likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);')
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);
disp(' ')

%Generate dataset input and underlying function parameters
disp('n = 20;')
xprev = linspace(0,4,201);
%r = [155,74,86,143,47,174,101,90,79,77,58,201,62,42,36,64,158,176,109,54,144,102,183,125,181,45,1,34,190,126,71,99,123,127,108,27,75,46,55,171,24,152,96,37,114,51,17,172,29,20,5,13,118,95,19,18,136,76,138,200,11,66,16,39,82,182,84,65,141,31,164,189,33,135,10,56,185,59,7,97,88,191,89,192,188,3,9,173,30,40,117,21,100,177,157,167,32,180,159,130,153,199,170,119,81,8,187,129,193,134,53,140,38,80,69,197,128,168,132,73,178,92,156,87,4,48,83,41,68,179,142,23,161,112,107,60,186,169,6,133,120,35,149,78,91,195,151,124,49,103,98,113,184,70,162,111,165,175,50,93,12,150,110,154,139,63,94,115,194,131,104,25,22,61,147,160,196,121,28,198,44,166,57,72,43,146,105,106,67,148,85,145,15,2,163,122,14,52,137,116,26];
%r = [14     5    17     3    10     4     8    19     7    21     6     9     1    13     2    16    20    15    11    12    18];
r = randperm(length(xprev));
x = xprev(r);
n = length(xprev);
l = floor(n*1/3);                   %size of training set
freq = 12;                          %frequency for underlying function
amp = .75;                            %amplitude for underlying function

traininds = [1:l];                  %indicies for different parts of data
testinds = [l+1:n];
testinds = [traininds testinds]; %to test over the whole distribution

%Define the underlying function and output w/ noise
%disp('x = gpml_randn(0.3, n, 1);')
%x = gpml_randn(0.3, n, 1);
disp('K = feval(covfunc{:}, hyp.cov, x);')
K = feval(covfunc{:}, hyp.cov, x(traininds)');
disp('mu = feval(meanfunc{:}, hyp.mean, x);')
mu = feval(meanfunc{:}, hyp.mean, x(traininds)');
disp('ytr = chol(K)''*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);')
%ytr = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);
ftr = amp*sin(freq*x(traininds).^.55)';     %training function - real underlying
%ftr = (1/3).*traininds';                    %training function - real underlying
noise = randn(length(x(traininds)),1);      %gaussian noise
ytr = ftr + .5*noise;                        %training output w/ noise

%Plot dataset and underlying function
figure(1)
set(gca, 'FontSize', 24)
disp(' '); disp('plot(x(traininds), ytr, ''+'')')
hold on;
plot(x(traininds), ytr, '+', 'MarkerSize', 12)
plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
axis([xprev(1) xprev(end) -2 2])
grid on
xlabel('input, x')
ylabel('output, ytr')
if write_fig, print -depsc f1.eps; end
disp(' '); disp('Hit any key to continue...'); pause

%Generate gaussian process w/ training data
disp(' ')
disp('nlml = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr)')
nlml = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds)', ytr)
disp(' ')

%Generate prediction w/ training data and test input
disp('z = linspace(-1.9, 1.9, 101)'';')
%z = linspace(0, 4, 201)';
z = xprev';
disp('[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds)', ytr, z);

%Plot initial prediction
figure(2)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];') 
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8);')
fill([z; flipdim(z,1)], f, [7 7 7]/8);

disp('hold on; plot(z, m); plot(x, ytr, ''+'')')
hold on; plot(z, m, 'LineWidth', 2); plot(x(traininds)', ytr, '+', 'MarkerSize', 12)
plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
axis([xprev(1) xprev(end) -2 2])
grid on
xlabel('input, x')
ylabel('output, ytr')
if write_fig, print -depsc f2.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp(' ')
disp('covfunc = @covSEiso; hyp2.cov = [0; 0]; hyp2.lik = log(0.1);')
covfunc = @covSEiso; hyp2.cov = [0; 0]; hyp2.lik = log(0.1);
disp('hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x, ytr)')
hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x, ytr);
disp(' ')

disp('exp(hyp2.lik)')
exp(hyp2.lik)
disp('nlml2 = gp(hyp2, @infExact, [], covfunc, likfunc, x, ytr)')
nlml2 = gp(hyp2, @infExact, [], covfunc, likfunc, x, ytr)
disp('[m s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x(traininds)', ytr, z);

%Optimizing the hyperparameters
disp(' ')
figure(3)
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)');
fill([z; flipdim(z,1)], f, [7 7 7]/8)
disp('hold on; plot(z, m); plot(x, ytr, ''+'')');
hold on; plot(z, m, 'LineWidth', 2); plot(x(traininds)', ytr, '+', 'MarkerSize', 12)
plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('input, x')
ylabel('output, ytr')
axis([xprev(1) xprev(end) -2 2])
if write_fig, print -depsc f3.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp(' ')
disp('hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);')
hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);
disp('hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, ytr);')
hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x(traininds)', ytr);
disp('[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds)', ytr, z);

%Refine the prediction
figure(4)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];')
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)')
fill([z; flipdim(z,1)], f, [7 7 7]/8)
disp('hold on; plot(z, m); plot(x, ytr, ''+'');')
hold on; plot(z, m, 'LineWidth', 2); plot(x(traininds)', ytr, '+', 'MarkerSize', 12)
plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('input, x')
ylabel('output, ytr')
axis([xprev(1) xprev(end) -2 2])
if write_fig, print -depsc f4.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp('large scale regression using the FITC approximation')
disp('nu = fix(n/2); u = linspace(.6,3.4,nu)'';')
nu = fix(n/2); u = linspace(.6,3.4,nu)';
disp('covfuncF = {@covFITC, {covfunc}, u};')
covfuncF = {@covFITC, {covfunc}, u};
disp('[mF s2F] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x, ytr, z);')
[mF s2F] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x(traininds)', ytr, z);

%Adjust training bounds
figure(5)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [mF+2*sqrt(s2F); flipdim(mF-2*sqrt(s2F),1)];')
f = [mF+2*sqrt(s2F); flipdim(mF-2*sqrt(s2F),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)')
fill([z; flipdim(z,1)], f, [7 7 7]/8)
disp('hold on; plot(z, mF); plot(x, ytr, ''+'');')
hold on; plot(z, mF, 'LineWidth', 2); plot(x(traininds)', ytr, '+', 'MarkerSize', 12)
disp('plot(u,1,''o'')')
plot(u,1,'ko', 'MarkerSize', 12)
plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('input, x')
ylabel('output, ytr')
axis([xprev(1) xprev(end) -2 2])
if write_fig, print -depsc f5.eps; end
disp(' ')