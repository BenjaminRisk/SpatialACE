function [results] = mle_ACE(y, X, familyst, initvalues)
% Benjamin Risk
% 6 March 2017
% Calculate MLE for a single location in the FSEM formulation of ACE model.
%
%INPUT:
% y: N x 1 vector
% X: N x p vector
% familyst: cell with
    % MZtp1 -- N x 1: ==1 if subject is the first twin for MZ pair
    % MZtp2 -- N x 1: ==1 if subject is the second twin for MZ pair
    % DZtp1 -- N x 1: ==1 if subject is the first twin for DZ pair
    % DZtp2 -- N x 1: ==1 if subject is the second twin for DZ pair
    % MDti -- N x 1: ==1 if subject is a singleton.
% initvalues -- (p+3) x 1 : optional initial values.
%OUTPUT:
%results structure containing:
% beta: vector of length p
% sigmasqA
% sigmasqC 
% sigmasqE
% yhat: N x 1 vector of predicted values
% resid: N x 1 vector of residuals
% f: -2*loglik at the estimated MLE
% exitflag_a: Optimizer output from alternative hypothesis 
%   1 Magnitude of gradient smaller than the TolFun tolerance.
%   2 Change in x was smaller than the TolX tolerance.
%   3 Change in the objective function value was less than the TolFun tolerance.
%   5 Predicted decrease in the objective function was less than the TolFun tolerance.
%   0 Number of iterations exceeded MaxIter or number of function evaluations exceeded MaxFunEvals.
%   -1 Algorithm was terminated by the output function.
%   -3 Objective function at current iteration went below ObjectiveLimit.
% lrt: likelihood ratio statistic for null: SigmaA=0 (-2*loglik(null) -
% (-2*loglik(alternative))
% pvalue: p-value of lrt
% exitflag_n: optimizer output from null hypothesis


MZtp1 = familyst.MZtp1;
MZtp2 = familyst.MZtp2;
DZtp1 = familyst.DZtp1;
DZtp2 = familyst.DZtp2;
MDti = familyst.MDti;


options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton','MaxIter',10000,'MaxFunEvals',1000*length(initvalues));
%options = optimset('GradObj','off', 'Display','off','Algorithm','quasi-newton');

[paramsA,f,exitflag] = fminunc(@(params) loglik_ACE(params, y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti), ...
           initvalues, options);
        
% [paramsA,f,exitflag] = fminunc(@(params) loglik_ACE(params, y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti), ...
%             initvalues);
        
results.betas = paramsA(4:end);
sigmasq = exp(paramsA(1:3));
results.sigmasqA = sigmasq(1);
results.sigmasqC = sigmasq(2);
results.sigmasqE = sigmasq(3);
yhattemp = X*results.betas;
results.fixefresid = y - yhattemp;
results.f = f;
results.exitflag_a = exitflag;

results2 = mle_CE(y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti, initvalues(2:end));
lrt = results2.f - f;
if lrt<0 
    lrt=0;
end
results.lrt = lrt;
results.pvalue = 1 - normcdf(sqrt(lrt));
results.exitflag_n = results2.exitflag;

% Calculate BLUPs:
N = size(X,1);
results.ai = zeros(N,1);
results.aij = zeros(N,1);
results.ci = zeros(N,1);

% MZs:
n1 = sum(MZtp1);
tindex = logical(MZtp1+MZtp2);
tempSigma = kron(eye(n1),inv(results.sigmasqA*ones(2)+results.sigmasqC*ones(2)+results.sigmasqE*eye(2)));
interblups = tempSigma*results.fixefresid(tindex);
results.ai(tindex) = kron(eye(n1),results.sigmasqA*ones(2))*interblups;
results.ci(tindex) = kron(eye(n1),results.sigmasqC*ones(2))*interblups;

% DZs:
n2 = sum(DZtp1);
tindex = logical(DZtp1+DZtp2);
tempSigma = kron(eye(n2),inv(0.5*results.sigmasqA*eye(2)+0.5*results.sigmasqA*ones(2)+...
    results.sigmasqC*ones(2)+results.sigmasqE*eye(2)));
interblups = tempSigma*results.fixefresid(tindex);
results.aij(tindex) = kron(eye(n2),0.5*results.sigmasqA*eye(2))*interblups;
results.ai(tindex) = kron(eye(n2),0.5*results.sigmasqA*ones(2))*interblups;
results.ci(tindex) = kron(eye(n2),results.sigmasqC*ones(2))*interblups;

% Singletons:
%n3 = sum(MDti);
%tempSigma = eye(n3)./(results.sigmasqA+results.sigmasqC+results.sigmasqE);
interblups = results.fixefresid(MDti)./(results.sigmasqA+results.sigmasqC+results.sigmasqE);
results.ai(MDti) = results.sigmasqA*interblups;
results.ci(MDti) = results.sigmasqC*interblups;

results.yhat = yhattemp+results.ai+results.aij+results.ci;
results.yfixed = yhattemp;
results.mixefresid = y - results.yhat;
% calculate approximate p-values for covariates:

% confidence bands for beta

% note: loglik_ACE defines -2*loglik

    tmp = hessian(@(params) loglik_ACE(params,y,X, familyst.MZtp1, familyst.MZtp2, familyst.DZtp1,...
            familyst.DZtp2, familyst.MDti),paramsA); %/(n1+n2+sum(MDti));


    %paramsA: (3+p) x 1 vector
    % WARNING: Shikai's code was incorrect
   warning('off','MATLAB:singularMatrix');
   warning('off','MATLAB:nearlySingularMatrix');
    tmp = diag(inv(tmp)); 
   warning('on','MATLAB:singularMatrix');
   warning('on','MATLAB:nearlySingularMatrix');
    % apply delta method to variance parameters
    % this version is correct, as determined from examining the
    % distribution of pvalues under the null hypothesis
    paramsA_se = sqrt(2)*sqrt(tmp).*[exp(paramsA(1:3));ones(length(paramsA)-3,1)];
    % (don't end up using the se's of variance components but retain for
    % future)
    % note the factor of 2 corrects for the fact that the likelihood
    % function is actually for -2*loglik
    results.betas_SE = paramsA_se(4:end);
    
    % TO DO: Replace with proportions correcting for finite sample size
    results.betas_pvalues = 2*(1-normcdf(abs(results.betas./results.betas_SE)));
    
%     % WARNING: I now only invert the fixed effects because variance
%     % components are sometime near zero leading to numeric instability.
%     tmp = diag(inv(tmp(4:end,4:end))); 
%     paramsA_se = 2*sqrt(tmp);
%     
 
%     % note the factor of 2 corrects for the fact that the likelihood
%     % function is actually for -2*loglik
%     results.betas_SE = paramsA_se;
%     results.betas_pvalues = 2*(1-normcdf(abs(results.betas./results.betas_SE)));
%     %results.betas_pvalues = normcdf(results.betas./results.betas_SE);
    
end

