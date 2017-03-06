function [results] = mle_CE(y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti,initvalues)
% Benjamin Risk
% Calculate MLE for a single location in the FSEM formulation of CE model,
% i.e., under H0: \sigma_a^2 = 0.
%
%Input:
% y: N x 1 vector
% X: N x p vector
% MZtp1 -- N x 1: ==1 if subject is the first twin for MZ pair
% MZtp2 -- N x 1: ==1 if subject is the second twin for MZ pair
% DZtp1 -- N x 1: ==1 if subject is the first twin for DZ pair
% DZtp2 -- N x 1: ==1 if subject is the second twin for DZ pair
% MDti -- N x 1: ==1 if subject is a singleton.
% initvalues: (2+p): optional input of initial values. 
%
%Output:
%results structure containing:
% betas: vector of length p
% sigmasq: vector of length 2: sigmasq_c, sigmasq_e
% yhat: N x 1 vector of predicted values
% resid: residuals, N x 1
% f: -2*loglik at the estimated MLE
% exitflag: 
%   1 Magnitude of gradient smaller than the TolFun tolerance.
%   2 Change in x was smaller than the TolX tolerance.
%   3 Change in the objective function value was less than the TolFun tolerance.
%   5 Predicted decrease in the objective function was less than the TolFun tolerance.
%   0 Number of iterations exceeded MaxIter or number of function evaluations exceeded MaxFunEvals.
%   -1 Algorithm was terminated by the output function.
%   -3 Objective function at current iteration went below ObjectiveLimit.

% if nargin<8 
%     p = size(X,2);
%     initvalues = zeros(p+2,1);
% end

options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton','MaxIter',10000,'MaxFunEvals',1000*length(initvalues));
[paramsN,f,exitflag] = fminunc(@(params) loglik_CE(params, y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti), ...
            initvalues, options);
%[paramsN,f,exitflag] = fminunc(@(params) loglik_CE(params, y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti), ...
%            initvalues);

results.betas = paramsN(3:end);
results.sigmasq = exp(paramsN(1:2));
results.yhat = X*results.betas;
results.resid = y - results.yhat;
results.f = f;
results.exitflag = exitflag;
end
