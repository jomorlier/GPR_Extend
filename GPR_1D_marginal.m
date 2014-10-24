function [f C] = GPR_1D_marginal()
% %Choose kernel
% kern = 1
% switch kern
%     case 1
% end
k = @(x,y) exp(-100*(x-y)'*(x-y));

addpath('..');

%Choose sample
%x = linspace(0,10,200);
%r = randperm(length(x));
%x = x(r);
x = [0:.05:.995];                  %x-input
r = [6 3 16 11 7 17 14 8 5 19 15 1 2 4 18 13 9 20 10 12];
x = x(r);
n = length(x);                      %size of total
l = floor(n*2/3);                   %size of training set
H = pseudo_kern_trick(x');
%Construct covariance matrix
C = zeros(n,n);
for i = 1:n
    for j = 1:n
        C(i,j) = k(x(i),x(j));
    end
end
ytr = mvnrnd(zeros(l,1),C(1:l,1:l)+.1*eye(l))';
beta = zeros(size(H));

%Construct covariance matrix blocks
% for i = 1:l
%     for j = 1:l
%         Kxx(i,j) = C(i,j);
%     end
% end
% for i = l+1:n
%     for j = 1:l
%         Kxox(i-l,j) = C(i,j);
%     end
% end
% for i = 1:l
%     for j = l+1:n
%         Kxxo(i,j-l) = C(i,j);
%     end
% end
% for i = l+1:n
%     for j = l+1:n
%         Kxoxo(i-l,j-l) = C(i,j);
%     end
% end
Kxx(1:l,1:l) = C(1:l,1:l) + .1*eye(length(C(1:l,1:l)));
Kxox(1:n-l,1:l) = C(l+1:end,1:l);
Kxxo(1:l,1:n-l) = C(1:l,l+1:end);
Kxoxo(1:n-l,1:n-l) = C(l+1:end,l+1:end);
%ftr = sin(x(1:l)/10)';

%Joint Posterior Distribution Computations
mu_all = C(:,1:l)*inv(Kxx)*ytr;
mu = Kxox*inv(Kxx)*ytr;         %New mean
D = Kxoxo - Kxox*inv(Kxx)*Kxxo; %New Covariance
D_all = C - C(:,1:l)*inv(Kxx)*C(1:l,:);

f = mvnrnd(mu,D,1);
%Sample from Gaussian Process at points
%u = randn(n-l,1);               %sample from normal disctribution
%[A,S,B] = svd(D);               %factor C
%f = A*sqrt(S)*u;                %resulting output

%Plot
%figure;
hold all;
clf;
plot(x(1:l),ytr,'o');
plot(x(l+1:end),mu,'c.');
plot(x,mu_all,'m.');
plot(x(l+1:end),f,'k.');

plot(repmat(x(l+1:end)',[1 size(f,1)])',f,'.');
axis([0,1,-2,2]);
end