function [rankest] = estrank(outAorC,prop)
%estrank estimate the number of components to include in the
%contrained covariance estimator
%INPUT:
%   outAorC: diagonal matrix of eigenvalues output from sandwich covariance
%   estimator for additive genetic or common environmental components
%   prop: proportion of variance to retain, e.g., 0.95

% one = diag(outAorC);
  one = outAorC;
  one(outAorC<0) = 0;
  two = cumsum(one);
  three = two./two(end);
  rankest = sum(three<prop);
end

