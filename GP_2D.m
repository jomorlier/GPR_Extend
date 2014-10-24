% %Choose kernel
% kern = 1
% switch kern
%     case 1
% end
k = @(x,y) exp(-100*(x-y)'*(x-y));

%Choose sample
points = (0:.05:1)';
[U,V] = meshgrid(points,points);
x = [U(:) V(:)]';
n = size(x,2);

%Construct covariance matrix
C = zeros(n,n);
for i = 1:n
    for j = 1:n
        C(i,j) = k(x(i),x(j));
    end
end

%Sample from Gaussian Process at points
u = randn(n,1);     %sample from normal distribution
[A,S,B] = svd(C);   %factor C
f = A*sqrt(S)*u;    %resulting output

%Plot
figure(2);
clf;
Z = reshape(f,sqrt(n),sqrt(n));
surf(U,V,Z);