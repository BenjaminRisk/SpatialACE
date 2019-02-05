function [f] = wtdloglik_CE(log_sigmaN, sumSR, sumXR, sumSMD, sumwts, n1, n2, n3)

% Loosely based on BLINDED FSEM_wlln 
% Calculate wtd loglik of variance components given residuals for CE model
% Most of my changes are efficiency improvements to prevent recalculating
% wtd residual sums for the same set of distances, etc. This function is designed to be
% used with an optimizer, as in mwle_CE
%Input:
% log_sigmaN: 2 x 1 vector
% sumSR see mwle_CE
% sumXR
% sumSMD
% n1: number of pairs of MZs
% n2: number of pairs of DZs
% n3: number of singletons
%
%Output:
% f: -2 weighted log likelihood under ACE model

% transformation so that sigma's are non-negative
sigma_c2 = exp(log_sigmaN(1));
sigma_e2 = exp(log_sigmaN(2));

a = sigma_c2+sigma_e2; 
b = sigma_c2;

f = ( (a*sumSR-b*sumXR)/(a^2-b^2) + (n1+n2)*log(a^2-b^2)*sumwts + ...
    sumSMD/a + n3*log(a)*sumwts )/sumwts;
end
