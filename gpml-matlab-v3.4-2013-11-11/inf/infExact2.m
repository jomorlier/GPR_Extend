function [post nlZ dnlZ] = infExact2(hyp, mean, cov, lik, x, y)

% Exact inference for a GP with Gaussian likelihood. Compute a parametrization
% of the posterior, the negative log marginal likelihood and its derivatives
% w.r.t. the hyperparameters. See also "help infMethods".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-06-10.
%                                      File automatically generated using noweb.
%
% See also INFMETHODS.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood');
end

[n, D] = size(x);


K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss
if sn2<1e-6                        % very tiny sn2 can lead to numerical trouble
  L = chol(K+sn2*eye(n)); sl =   1;   % Cholesky factor of covariance with noise
  pL = -solve_chol(L,eye(n));                            % L = -inv(K+inv(sW^2))
else
  L = chol(K/sn2+eye(n)); sl = sn2;                       % Cholesky factor of B
  pL = L;                                           % L = chol(eye(n)+sW*sW'.*K)
end

    if xeqz
        q = x;
    else
        q = [x; z];
    end
    H = [sin(x) cos(x)]';

    Hq = [sin(q) cos(q)]';
    
    b = [-1; -1];
    sf2 = exp(2*hyp(end));
    Ky = Ko+sf2*eye(size(Ko));
    sf2 = exp(2*hyp(end)); 
    sigma_B = hyp(end-1); % possible that this should be the log to follow the other hyppars
    N = size(H,1);
    B = sigma_B^2*eye(N);
    B_inv = inv(B);
    deriv_B = 2*sigma_B*eye(N); 
    Ky_inv = inv(Ky);
    A = inv(B) + H*Ky_inv*H';
    A_inv = inv(A);
    
    R = Hq - H*Ky_inv*Kq;
    K = Kqq + R'*inv(B_inv+H*Ky_inv*H')*R;  %with xeqz, symmetric matrix kxx
    
    if ~xeqz %cross covariances kxz
        K = K(length(x)+1:end,1:length(x))';
    end
    if nargout == 2
        beta_mean = inv(B_inv+H*Ky_inv*H');
        beta_mean = beta_mean*(H*Ky_inv*y+B_inv*b);
        Basis_Addend = R'*beta_mean;
        Basis_Addend = Basis_Addend(length(x)+1:end);
    end
    
    
alpha = solve_chol(L,y-m)/sl;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                  % sqrt of noise precision vector
post.L = pL;

if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sl)/2;   % -log marg lik
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;                                 % allocate space for derivatives
    Q = solve_chol(L,eye(n))/sl - alpha*alpha';     % precompute for convenience
    for i = 1:numel(hyp.cov)
        dnlZ.cov(i) = sum(sum(Q.*feval(cov{:}, hyp.cov, x, [], i)))/2;
    end
    dnlZ.lik = sn2*trace(Q);
    for i = 1:numel(hyp.mean), 
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*alpha;
    end
  elseif i==length(v)+2
      K = -1/2*y'*Ky_inv*H'*A_inv*B_inv*deriv_B*B_inv*A_inv*H*Ky_inv*y -1/2*trace(B_inv*deriv_B) -1/2*trace(-A_inv*B_inv*deriv_B*B_inv);
  end
  end
end
