function [results] = mwle_ACE(R, dist, h, familyst, initvalues)
% Calculate weighted MLE for a single location in the FSEM formulation of ACE model.
%
%Input:
% R:    N x M matrix of residuals (M locations)
% dist: M x 1 vector of distances from focal location
% h:    bandwidth for the biweight kernel; note that all vertices with
%       dist>=h have weight=0
% familyst
%   MZtp1 -- N x 1: ==1 if subject is the first twin for MZ pair
%   MZtp2 -- N x 1: ==1 if subject is the second twin for MZ pair
%   DZtp1 -- N x 1: ==1 if subject is the first twin for DZ pair
%   DZtp2 -- N x 1: ==1 if subject is the second twin for DZ pair
% MDti -- N x 1: ==1 if subject is a singleton.
% initvalues -- 3 x 1 : optional initial values. on log scale.
%Output:
% sigmasq: vector of length 3; sigmasq_a, sigmasq_c, sigmasq_e
% f: -2*loglik at the estimated MWLE
% exitflag: 
%   1 Magnitude of gradient smaller than the TolFun tolerance.
%   2 Change in x was smaller than the TolX tolerance.
%   3 Change in the objective function value was less than the TolFun tolerance.
%   5 Predicted decrease in the objective function was less than the TolFun tolerance.
%   0 Number of iterations exceeded MaxIter or number of function evaluations exceeded MaxFunEvals.
%   -1 Algorithm was terminated by the output function.
%   -3 Objective function at current iteration went below ObjectiveLimit.

if nargin<5
    initvalues = zeros(3,1);
end

MZtp1 = familyst.MZtp1;
MZtp2 = familyst.MZtp2;
DZtp1 = familyst.DZtp1;
DZtp2 = familyst.DZtp2;
MDti = familyst.MDti;

[nSubject,nVertex] = size(R);
if nSubject>nVertex
    warning('R has more rows than columns -- check that it is N x M')
end
index = (dist<=h);
subR = R(:,index);
wts = biweight(dist(index)/h)/h;
sumwts = sum(wts);

% Wtd residuals do not change so calculate outside of optimizer:
n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);
SM = subR(MZtp1,:).^2 + subR(MZtp2,:).^2; XM = 2*subR(MZtp1,:).*subR(MZtp2,:); % n1 x M
SD = subR(DZtp1,:).^2 + subR(DZtp2,:).^2; XD = 2*subR(DZtp1,:).*subR(DZtp2,:); % n2 x M
SMD = subR(MDti,:).^2; % n3 x M

sumSM = sum(SM)*wts;
sumXM = sum(XM)*wts;
sumSD = sum(SD)*wts;
sumXD = sum(XD)*wts;
sumSMD = sum(SMD)*wts;

%options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton');
options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton','MaxIter',10000,'MaxFunEvals',1000*length(initvalues));
[paramsA,f,exitflag] = fminunc(@(params) wtdloglik_ACE(params, sumSM, sumXM, sumSD, sumXD, sumSMD, sumwts, n1, n2, n3), ...
            initvalues, options);
sigmasq = exp(paramsA);
results.sigmasqA = sigmasq(1);
results.sigmasqC = sigmasq(2);
results.sigmasqE = sigmasq(3);
results.f = f;
results.exitflag_a = exitflag;

end
