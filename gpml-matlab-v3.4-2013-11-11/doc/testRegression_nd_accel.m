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

[sweep,sweep_vecs] = generate_sweep; %Generate some test input

%Generate dataset input and underlying function parameters
[sim_ds,E] = orig_sim;%orig_sim(0.01,100,pi/12,1,1.05,0,pi/2,0);
xprev = sim_ds(1:(end/30),2:end);
timeset = sim_ds(1:(end/30),1);
%r = randperm(length(xprev(:,1)));
r = [798,499,470,356,173,920,819,524,881,183,72,472,163,581,139,771,792,825,193,220,233,545,28,88,16,123,549,234,941,676,536,351,948,580,256,109,405,703,122,993,753,279,152,490,712,829,672,164,782,959,612,926,413,589,517,219,459,811,879,331,725,976,613,947,774,388,762,738,913,915,272,853,627,943,58,812,815,17,140,515,240,227,795,567,93,14,347,375,818,112,50,733,352,374,875,340,719,31,38,610,145,531,494,354,638,551,159,648,917,970,384,523,414,294,378,925,802,585,430,604,20,420,237,749,75,860,527,911,950,29,932,453,360,496,626,298,51,635,3,396,601,546,509,128,7,23,843,306,26,824,379,464,906,885,320,880,974,513,641,902,452,170,497,254,983,160,852,493,969,615,845,283,19,519,172,653,599,313,87,155,495,41,245,748,867,471,556,84,377,325,968,66,886,758,955,571,9,156,399,768,268,526,304,359,877,104,432,131,579,639,415,557,645,296,995,942,369,597,598,862,940,308,117,781,199,688,324,530,376,502,44,731,591,54,62,582,696,870,315,48,813,100,723,864,76,49,357,540,425,252,390,205,249,207,446,80,63,456,1000,134,221,419,346,403,904,541,265,742,734,166,45,192,68,264,345,142,141,255,434,189,892,394,732,305,846,789,36,422,327,756,514,586,480,146,412,278,113,368,466,806,922,35,997,662,74,713,784,132,717,695,323,218,65,105,481,350,763,785,126,410,873,246,674,211,842,448,777,794,757,372,309,737,779,223,963,625,923,992,692,791,752,525,355,373,685,103,889,501,401,722,232,735,958,721,444,741,898,831,949,276,167,887,871,708,553,408,176,342,851,235,698,655,120,32,348,133,1,869,230,764,999,914,150,455,966,839,709,57,930,267,686,433,292,11,951,584,2,787,882,143,191,182,847,321,96,251,972,6,437,799,289,380,578,196,537,624,431,986,977,94,664,416,442,83,683,259,338,908,522,4,894,478,707,81,106,406,689,258,27,482,652,903,349,876,79,288,195,854,34,554,814,521,823,552,367,286,91,729,603,936,460,344,206,975,445,594,705,661,391,632,840,590,400,370,761,418,398,441,409,504,693,766,13,697,301,617,184,365,263,668,747,861,726,290,447,772,198,37,775,952,874,808,572,144,197,602,780,836,678,439,435,78,248,790,200,275,277,665,884,595,728,716,618,202,242,397,171,918,224,70,67,651,161,175,426,60,225,675,467,859,244,834,901,71,991,138,534,99,312,803,82,273,209,820,727,226,535,856,449,228,392,793,921,934,454,751,210,488,569,767,560,311,341,307,487,440,18,817,699,844,208,314,508,169,649,153,924,516,30,335,663,865,770,25,485,620,465,361,804,411,451,701,529,670,611,168,260,667,609,339,300,985,121,486,575,945,555,633,387,257,424,755,111,130,570,654,47,15,896,190,95,743,343,56,187,124,588,563,642,532,691,897,125,681,680,280,596,788,899,217,559,544,533,826,450,989,457,383,891,905,318,363,800,659,711,101,765,148,90,366,646,682,158,330,270,883,498,266,262,402,404,186,636,85,382,973,677,539,61,998,135,739,274,606,643,573,987,423,46,631,436,931,303,927,754,194,241,907,608,468,730,505,619,592,503,243,282,151,759,69,185,107,229,954,463,714,53,473,961,801,238,214,933,212,371,710,706,900,838,518,607,939,511,236,909,86,750,916,461,850,855,326,938,868,203,562,295,605,796,462,833,333,215,805,890,715,816,872,893,981,979,957,773,679,669,700,746,97,293,666,538,328,22,644,937,629,965,291,960,512,500,322,622,736,59,489,944,395,744,593,623,284,386,180,616,458,253,994,429,364,385,77,953,261,797,786,353,657,807,720,154,73,809,574,428,492,769,40,10,830,421,358,630,157,477,64,336,407,332,165,821,858,92,317,702,984,174,718,362,929,740,690,102,547,127,507,149,848,39,565,849,583,137,684,928,98,89,381,188,239,484,33,329,543,978,660,946,895,982,443,962,866,990,971,55,566,671,129,778,760,285,201,837,281,510,996,319,119,506,658,564,827,956,561,637,110,389,299,491,302,647,271,568,427,43,863,483,832,980,910,878,967,835,474,520,115,810,269,776,108,438,287,783,745,600,888,704,964,136,179,634,213,5,828,310,628,479,316,857,724,337,204,841,222,114,8,614,42,118,21,935,181,24,919,231,822,216,673,250,147,475,912,162,550,988,469,334,116,548,393,247,640,542,694,576,656,650,417,52,177,297,12,558,621,687,577,476,528,178,587];
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
ftr = diff(xprev(:,2))';                         %outputs are derivatives
ftr = ftr(1,r(traininds));
noise = randn(1,length(ftr));                      %gaussian noise
ytr = ftr' + .005*noise';                         %training output w/ noise

%Plot dataset and underlying function
figure(2)
set(gca, 'FontSize', 24)
disp(' '); disp('plot(x(traininds), ytr, ''+'')')
hold on;
plot(x(traininds,1), ytr(traininds,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
axis([min(x(:,1))-2 max(x(:,1))+2 -.3 .3])
grid on
xlabel('input x')
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
%z = sweep;
z = [linspace(-2*pi,2*pi,1000)' linspace(-4,4,1000)' linspace(-5,5,1000)'];
disp('[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, ytr, z);')
[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x(traininds,:), ytr, z);

%Plot initial prediction
figure(3)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];') 
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8);')
%fill([z(:,1); flipdim(z(:,1),1)], f, [7 7 7]/8)

disp('hold on; plot(z, m); plot(x, ytr, ''+'')')
hold on; plot(z(:,1), '.', 'Markersize',12);
plot(x(traininds,1), ytr(traininds,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
axis([min(x(:,1))-2 max(x(:,1))+2 -.3 .3])
grid on
xlabel('input x')
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
%fill([z(:,1); flipdim(z(:,1),1)], f, [7 7 7]/8)
disp('hold on; plot(z, m); plot(x, ytr, ''+'')');
hold on; plot(z(:,1), m, '.', 'Markersize',12);

plot(x(traininds,1), ytr(traininds,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('input x')
ylabel('output, ytr')
axis([min(x(:,1))-2 max(x(:,1))+2 -.3 .3])
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
%fill([z(:,1); flipdim(z(:,1),1)], f, [7 7 7]/8)
disp('hold on; plot(z, m); plot(x, ytr, ''+'');')
hold on; plot(z(:,1), m, '.', 'Markersize',12);

plot(x(traininds,1), ytr(traininds,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('input x')
ylabel('output, ytr')
axis([min(x(:,1))-2 max(x(:,1))+2 -.3 .3])
if write_fig, print -depsc f4.eps; end
disp(' ');

disp('large scale regression using the FITC approximation')
disp('nu = fix(n/2); u = linspace(.6,3.4,nu)'';')
nu = fix(n/2); pu = linspace(.6,3.4,nu)'; vu = linspace(.6,3.4,nu)'; uu = linspace(.6,3.4,nu)'; u = [pu,vu,uu]; 
disp('covfuncF = {@covFITC, {covfunc}, u};')
covfuncF = {@covFITC, {covfunc}, u(traininds,:)};
disp('[mF s2F] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x, ytr, z);')
[mF_accel s2F] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x(traininds,:), ytr, z);

%Adjust training bounds
figure(6)
set(gca, 'FontSize', 24)
disp(' ')
disp('f = [mF+2*sqrt(s2F); flipdim(mF-2*sqrt(s2F),1)];')
f = [mF_accel+2*sqrt(s2F); flipdim(mF_accel-2*sqrt(s2F),1)];
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8)')
%fill([z(:,1); flipdim(z(:,1),1)], f, [7 7 7]/8)
disp('hold on; plot(z, mF); plot(x, ytr, ''+'');')
hold on; plot(z(:,1), mF_accel, '.', 'Markersize',12);

plot(x(traininds,1), ytr(traininds,1), 'r+', 'MarkerSize', 12)
%plot(x(traininds,2), ytr(:,2), 'g+', 'MarkerSize', 12)
disp('plot(u,1,''o'')')
plot(pu,1,'ko', 'MarkerSize', 12)
%plot(xprev,amp*sin(freq*xprev'.^.55),'k-');
grid on
xlabel('input x')
ylabel('output, ytr')
axis([min(x(:,1))-2 max(x(:,1))+2 -.3 .3])
if write_fig, print -depsc f5.eps; end
disp(' ')