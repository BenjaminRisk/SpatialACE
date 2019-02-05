function [f] = wtdloglik_ACE(log_sigmaA, sumSM, sumXM, sumSD, sumXD, sumSMD, sumwts, n1, n2, n3)

% Loosely based on BLINDED FSEM_wlla 
% Calculate wtd loglik of variance components given residuals
% Most of my changes are efficiency improvements to prevent recalculating
% wtd residual sums for the same set of distances, etc. This function is designed to be
% used with an optimizer, as in mwle_ACE
%Input:
% log_sigmaA: 3 x 1 vector
% sumSM see mwle_ACE
% sumXM
% sumSD
% sumXD
% sumSMD
% sumwts
% n1: number of pairs of MZs
% n2: number of pairs of DZs
% n3: number of singletons
%
%Output:
% f: -2 weighted log likelihood under ACE model

% transformation so that sigma's are non-negative
sigma_a2 = exp(log_sigmaA(1));
sigma_c2 = exp(log_sigmaA(2));
sigma_e2 = exp(log_sigmaA(3));

a = sigma_a2 + sigma_c2 + sigma_e2;
b = sigma_a2 + sigma_c2;
c = 0.5*sigma_a2 + sigma_c2;

%BRisk: note: I rescaled so wts sum to one 
% I also no longer divide by M...
% I also return -2*loglik instead of -loglik

%      f = ( (a*sum(SM(:))-b*sum(XM(:)))/(a^2-b^2) + n1*log(a^2-b^2)*sw + ...
%          (a*sum(SD(:))-c*sum(XD(:)))/(a^2-c^2) + n2*log(a^2-c^2)*sw + ...
%          sum(SMD(:))/a + n3*log(a)*sw )/2/M; % - loglikelihood 

f = ( (a*sumSM-b*sumXM)/(a^2-b^2) + n1*log(a^2-b^2)*sumwts + ...
         (a*sumSD-c*sumXD)/(a^2-c^2) + n2*log(a^2-c^2)*sumwts + ...
         sumSMD/a + n3*log(a)*sumwts )/sumwts; % - 2*loglikelihood 

end
