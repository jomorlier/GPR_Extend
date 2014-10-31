function K = covConst_2(hyp, x, z, i)
% dbstack
% Intermediate development cov function, to prepare for progressively more
% complex basis functions.  Reference made to covConst.
% hyperparameter (the constant) parameterized as itself (?) - does it need
% to be log(c^2)? - rather ,log(c^2)
if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; 
dg = strcmp(z,'diag') && numel(z)>0;        % determine mode  
c2 = exp(hyp(1)*2);
ell = 1;

% precompute squared distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sq_dist(x'/ell);
  else                                                   % cross covariances Kxz
    K = sq_dist(x'/ell,z'/ell); %-prediction
  end
  K = ones(size(K))*c2;
end

if nargin>3                                                        % covariances
    if i==1
        K = 2*K;
    else                                                        % derivatives
        error('Unknown hyperparameter')
    end
end