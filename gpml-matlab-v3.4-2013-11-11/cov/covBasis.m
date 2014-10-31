function K = covBasis(hyp, x, z, i)

% covBasis - implementation of suggested fixed basis function covariance
% sum from GPML Ch. 2.7, p. 28.  To be called with cov sum.
% Hyperparameters represent B, the coefficient(s) of the basis functions
% (see Eq. 2.39-2.42)


if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0;
dg = strcmp(z,'diag') && numel(z)>0;        % determine mode  
%if the third arg (z) is not diag, but exists, then the mode is prediction
%because z is the test points (xeqz T, dg F)

H = [sin(x) cos(x)];
Hz = [sin(z) cos(z)];
sigma_B = hyp(1)^2; % possible that this should be the log to follow the other hyppars
B = sigma_B^2*eye(N);
B_inv = inv(B);
deriv_B = 2*sigma_B*eye(N); 
Ky_inv = inv(Ky);
A = inv(B) + H*Ky_inv*H';
A_inv = inv(A);
if xeqz
    y = x;
else
    y = z;
end


% precompute squared distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sq_dist(x'/ell);
  else                                                   % cross covariances Kxz
    K = sq_dist(x'/ell,z'/ell); %-prediction
  end
end

if nargin<4                                                        % covariances
  K = sf2*exp(-K/2);
else                                                               % derivatives
  if i==1
    K = -1/2*y'*Ky_inv*H'*A_inv*B_inv*deriv_B*B_inv*A_inv*H*Ky_inv*y -1/2*trace(B_inv*deriv_B) -1/2*trace(-A_inv*B_inv*deriv_B*B_inv);
  else
    error('Unknown hyperparameter')
  end
end