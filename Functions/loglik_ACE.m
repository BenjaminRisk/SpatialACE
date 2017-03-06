function [f] = loglik_ACE(paramsA, y, X, MZtp1, MZtp2, DZtp1, DZtp2, MDti)
% Benjamin Risk
% Adapted from Shikai Luo's FSEM_lla 
% Computes the -2loglikelihood of the functional SEM formulation of the ACE model.
%
%Input:
% paramsA: (3+p) x 1 vector
% y: N x 1 vector -- MUST be sorted by numeric zygosity 
%   (1 = MZ, 2 = DZ, 3 = singleton) then familyID, as created in 
%   a2CreateSubjectMatrix.m.
% X: N x p matrix
% MZtp1 -- N x 1: ==1 if subject is the first twin for MZ pair
% MZtp2 -- N x 1: ==1 if subject is the second twin for MZ pair
% DZtp1 -- N x 1: ==1 if subject is the first twin for DZ pair
% DZtp2 -- N x 1: ==1 if subject is the second twin for DZ pair
% MDti -- N x 1: ==1 if subject is a singleton.
%
%Output 
% f: -2loglikelihood

% transformation so that sigma's are non-negative
sigma_a2 = exp(paramsA(1)); 
sigma_c2 = exp(paramsA(2)); 
sigma_e2 = exp(paramsA(3));

beta = paramsA(4:end);

R = y - X*beta;

a = sigma_a2+sigma_c2+sigma_e2; 
b = sigma_a2+sigma_c2; 
c = 0.5*sigma_a2+sigma_c2;

n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);

SM = sum(R(MZtp1).^2 + R(MZtp2).^2); 
XM = 2*sum(R(MZtp1).*R(MZtp2)); 

SD = sum(R(DZtp1).^2 + R(DZtp2).^2);
XD = 2*sum(R(DZtp1).*R(DZtp2));

SMD = sum(R(MDti).^2);

f = (a*SM-b*XM)/(a^2-b^2) + n1*log(a^2-b^2) + ...
    (a*SD-c*XD)/(a^2-c^2) + n2*log(a^2-c^2) + SMD/a + n3*log(a);

end