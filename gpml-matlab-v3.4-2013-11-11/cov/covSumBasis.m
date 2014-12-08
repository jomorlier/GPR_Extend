function [K, Basis_Addend] = covSumBasis(cov, hyp, x, z, i)
% dbstack
% covSumBasis - another attempt at the implementation of additional basis
% function code, from GPML 2.7 Eq. 2.39-42
% MBocamazo 2014-Oct

if numel(cov)==0, error('We require at least one summand.'), end
for ii = 1:numel(cov)                        % iterate over covariance functions
  f = cov(ii); 
  if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  j(ii) = cellstr(feval(f{:}));                          % collect number hypers
end

if nargin<3                                        % report number of parameters
  K = char(j(1)); 
  for ii=2:length(cov)
      K = [K, '+', char(j(ii))];
  end
  K = [K, '+', '2']; %two hyppars are beta (end-1) and sigma (noise) (end)
  return
end

if nargin<4, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; 
dg = strcmp(z,'diag') && numel(z)>0;
%if the third arg (z) is not diag, but exists, then the mode is prediction
%because z is the test points (xeqz T, dg F)

[n,D] = size(x);

v = [];               % v vector indicates to which covariance parameters belong
for ii = 1:length(cov)
    v = [v repmat(ii, 1, eval(char(j(ii))))]; 
end

if nargin<=6                                                        % covariances
    if dg
        K = zeros(size(x,1),1); return
    end
    K = 0;     
    for ii = 1:length(cov)                      % iteration over summand functions
        f = cov(ii); 
        if iscell(f{:}), f = f{:}; end % expand cell array if necessary
        K = K + feval(f{:}, hyp(v==ii), x, x);              % accumulate covariances
    end
    if xeqz
        q = x;
    else
        q = [x; z];
    end
    Ko = K; %save Ko
    K = 0;
    for ii = 1:length(cov)                      % iteration over summand functions
        f = cov(ii); 
        if iscell(f{:}), f = f{:}; end % expand cell array if necessary
        K = K + feval(f{:}, hyp(v==ii), x, q);              % accumulate covariances
    end
    Kq = K;
    
    K = 0;
    for ii = 1:length(cov)                      % iteration over summand functions
        f = cov(ii); 
        if                                                                                                  iscell(f{:}), f = f{:}; end % expand cell array if necessary
        K = K + feval(f{:}, hyp(v==ii), q, q);              % accumulate covariances
    end
    Kqq = K;
    
end
%need cov calc'ed to comp the derivative

if (nargin==6) && (nargout == 1)                                   % derivatives  
  if i<=length(v)
    vi = v(i);                                       % which covariance function
    j = sum(v(1:i)==vi);                    % which parameter in that covariance
    f  = cov(vi);
    if iscell(f{:}), f = f{:}; end         % dereference cell array if necessary
    K = feval(f{:}, hyp(v==vi), x, z, j);                   % compute derivative
  elseif i==length(v)+1
      K = sf2*exp(-K/2).*K;
  else
    error('Unknown hyperparameter')
  end
end