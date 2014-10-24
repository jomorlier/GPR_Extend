function K = covZero(hyp, x, z, i)
% dbstack
% Intermediate development covariance function.  Returns all zero
% covariance.  Takes zero hyperparameters.  modified from covSEiso so all
% the args still line up nicely. Michael Bocamazo 2014-10-24

if nargin<2, K = '0'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; 
dg = strcmp(z,'diag') && numel(z)>0;        % determine mode  
ell = 1; sf2 = 1;

% precompute squared distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sq_dist(x'/ell);
  else                                                   % cross covariances Kxz
    K = sq_dist(x'/ell,z'/ell); %-prediction
  end
  K = zeros(size(K));
end

if nargin<4                                                        % covariances
  K = sf2*exp(-K/2);
  K = zeros(size(K));
else                                                               % derivatives
    K = [];
    %question: [] or 0?
%   if i==1
%     K = sf2*exp(-K/2).*K;
%   elseif i==2
%     K = 2*sf2*exp(-K/2);
  if i ~= 0
    error('Unknown hyperparameter')
  end
end