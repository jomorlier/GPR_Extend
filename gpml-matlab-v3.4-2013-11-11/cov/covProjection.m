function K = covProjection(cov, hyp, x, z, i)

% Implementation of projection by basis functions into higher dimensions as
% seen in Ch. 2.7 of GPML.  Should act as a wrapper for other covariance
% functions and take x as struct, with the second field as H.  Based
% on covSum, similar bookkeeping.
% Michael Bocamazo, 2014/08/11

if numel(cov)~=1, error('covProjection requires exactly one covFcn'), end 
% Question for generality: could a composite be the one fcn? Then the
% following would still work:
for ii = 1:numel(cov)                        % iterate over covariance functions
  f = cov(ii); 
  if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  j(ii) = cellstr(feval(f{:}));                          % collect number hypers
end
% Else in the one arg case j(1) = covFcn
if nargin<3                                        % report number of parameters
  K = char(j(1)); 
  for ii=2:length(cov) %and this shouldn't break on an empty matrix if original works
      K = [K, '+', char(j(ii))];
  end, return
end

%if x isstruct here too
if nargin<4, z = []; end                                   % make sure, z exists

[n,D] = size(x); 

v = [];               % v vector indicates to which covariance parameters belong
for ii = 1:length(cov)
    v = [v repmat(ii, 1, eval(char(j(ii))))]; 
end

if isstruct(x)
    B = hyp(end); %this is a somewhat general, could be better
    hyp_cov = hyp(1:end-1);
%     x = x.raw;
%     z = z.raw;
    if nargin<5                                                        % covariances
      K = 0; 
      if nargin==3, z.raw = []; end                                 % set default
      for ii = 1:length(cov)                      % iteration over summand functions
        f = cov(ii); 
        if iscell(f{:}), f = f{:}; end % expand cell array if necessary
        K = K + feval(f{:}, hyp(v==ii), x.raw, z.raw);              % accumulate covariances        
      end    
      % now use basis fcns as well
      % fevals just run already handle the diag/non diag cases
      % contributions
      if strcmp(z.raw,'diag')
          %handle signal variance here
          %K = K + ? or K = ...
          covProj_contribution = feval(cov{:}, hyp.cov, xs(id,:), 'diag');
      elseif numel(z.raw)>0 && ~strcmp(z.raw,'diag') %this case doesn't make sense
          % self-covariance  (?)
          covProj_contribution = x.proj'*B*x.proj; %?
      else
          % return predictive covariance, the actual object of interest
          Kxx = feval(cov{:}, hyp.cov, x.raw); %self covariance
          Kxxo = feval(cov{:}, hyp.cov, x.raw, z.raw); %cross-covariance
          R = z.proj - x.proj*inv(Kxx)*Kxxo;
          covProj_contribution = R'*inv(inv(B)+x.proj*inv(Kxx)*x.proj')*R;          
      end      
      K = K + covProj_contribution;
    else                                                               % derivatives
      if i<=length(v)
        vi = v(i);                                       % which covariance function
        j = sum(v(1:i)==vi);                    % which parameter in that covariance
        f  = cov(vi);
        if iscell(f{:}), f = f{:}; end         % dereference cell array if necessary
        K = feval(f{:}, hyp(v==vi), x.raw, z.raw, j);                   % compute derivative
      else
        % basis function hyp must be handled here or near here, don't
        % hard-code hypProj numel check
        error('Unknown hyperparameter')
      end
    end
    
else %if x is not struct, same as covSum
    if nargin<5                                                        % covariances
      K = 0; 
      if nargin==3, z = []; end                                 % set default
      for ii = 1:length(cov)                      % iteration over summand functions
        f = cov(ii); 
        if iscell(f{:}), f = f{:}; end % expand cell array if necessary
        K = K + feval(f{:}, hyp(v==ii), x, z);              % accumulate covariances
      end
    else                                                               % derivatives
      if i<=length(v)
        vi = v(i);                                       % which covariance function
        j = sum(v(1:i)==vi);                    % which parameter in that covariance
        f  = cov(vi);
        if iscell(f{:}), f = f{:}; end         % dereference cell array if necessary
        K = feval(f{:}, hyp(v==vi), x, z, j);                   % compute derivative
      else
        error('Unknown hyperparameter')
      end
    end
end
