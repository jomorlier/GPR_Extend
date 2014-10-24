% %Choose kernel
% kern = 1
% switch kern
%     case 1
% end
k = @(x,y) exp(-100*(x-y)'*(x-y));

%Choose sample
x = (0:.005:1);
n = length(x);

%Construct covariance matrix
C = zeros(n,n);
for i = 1:n
    for j = 1:n
        C(i,j) = k(x(i),x(j));
    end
end

%Sample from Gaussian Process at points
u = randn(n,1);     %sample from normal distribution
[A,S,B] = svd(C);   %factor C with diagonal S and unitary matrices A, B
f = A*sqrt(S)*u;    %resulting output

%Plot
hold all;
plot(x,f,'.','LineWidth',2);
axis([0,1,-2,2]);