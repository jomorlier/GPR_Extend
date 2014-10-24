% %Choose kernel
% kern = 1
% switch kern
%     case 1
% end
rng('default')

%Covariance function
k = @(x,y) exp(-100*(x-y)'*(x-y));

%Choose sample
%x = linspace(0,10,200);
%r = randperm(length(x));
%x = x(r);

x = [0:.05:2];                  %x-input
r = randperm(length(x));
x = x(r);
n = length(x);                      %size of total
l = floor(n*1/2);                   %size of training set
ftr = 3*sin(5*x(1:l))';                    %training function - real underlying
noise = randn(length(x(1:l)),1);      %gaussian noise
ytr = ftr + noise;                  %training output w/ noise

%Construct covariance matrix
C = zeros(n,n);
for i = 1:n
    for j = 1:n
        C(i,j) = k(x(i),x(j));
    end
end
disp(C)

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
Kxx(1:l,1:l) = C(1:l,1:l) + .05*eye(length(C(1:l,1:l))); %coefficient encodes prior noise belief (?)
Kxox(1:n-l,1:l) = C(l+1:end,1:l); %noise is not included in any other covariance submatrix
Kxxo(1:l,1:n-l) = C(1:l,l+1:end);
Kxoxo(1:n-l,1:n-l) = C(l+1:end,l+1:end);
%ftr = sin(x(1:l)/10)';

%Joint Posterior Distribution Computations
mu_all = C(:,1:l)*inv(Kxx)*ytr; %train covariance, inv self cov, train -> mean all
mu = Kxox*inv(Kxx)*ytr;         %New mean
D = Kxoxo - Kxox*inv(Kxx)*Kxxo; %New Covariance ... difference from gaussian processes book: missing transpose?
D_all = C - C(:,1:l)*inv(Kxx)*C(1:l,:); 

f = mvnrnd(mu,D,50);
%Sample from Gaussian Process at points
%u = randn(n-l,1);               %sample from normal distribution
%[A,S,B] = svd(D);               %factor C
%f = A*sqrt(S)*u;                %resulting output

%Plot
%figure;
hold all;
%clf;
plot(x(1:l),ytr,'o');
plot(x(1:l),ftr,'ro');
% plot(x(l+1:end),mu,'co');
plot(x,mu_all,'m.');
%Note: labels unsure
legend('train with noise','train without noise','output')
% plot(repmat(x(l+1:end)',[1 size(f,1)])',f,'.');
%axis([0,1,-2,2]);