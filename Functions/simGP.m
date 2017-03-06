function [xmat] = simGP(covmat, rank, n)
%simGP Simulate a Gaussian Process or degenerate multivariate normal
%       distribution.
% Author: Benjamin Risk, 9/15/2015
%Input:
%   covmat: M x M covariance matrix
%   rank: rank of the covariance matrix. This function
%       is efficient when rank<<M for known rank
%   n: number of Gaussian processes to simulate, e.g., number of subjects
%
%Output:
%   xmat: M x n matrix where the finite dimensional distribution of each
%       row is equal to covmat

% since rank << M, eigs results in huge time savings
[psi,lambda] = eigs(covmat,rank,'la');
z = randn(rank,n);

%xmat = bsxfun(@times,diag(lambda),z); %this does not seem to result in
%speed improvements

xmat = sqrt(lambda)*z;
xmat = psi*xmat;

end

