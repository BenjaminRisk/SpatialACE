function [hvec,mse] = cv_mwle_ACE(R, dist, familyID, MZtp1, MZtp2, DZtp1, DZtp2, MDti, initvalues)
%cv_mwle_ACE Use cross-validation to select the bandwidth in
% the weighted mle for the variance paramaters
%Input:
% R:    N x V matrix of residuals (V locations, N individuals)
% dist: V x 1 vector of distances from focal location
% familyID -- N x 1 : integer-valued vector with the family id for each
%   subject
% MZtp1 -- N x 1: ==1 if subject is the first twin for MZ pair
% MZtp2 -- N x 1: ==1 if subject is the second twin for MZ pair
% DZtp1 -- N x 1: ==1 if subject is the first twin for DZ pair
% DZtp2 -- N x 1: ==1 if subject is the second twin for DZ pair
% MDti -- N x 1: ==1 if subject is a singleton.
% initvalues -- 3 x 1 : optional initial values on log scale.
%Output:
% hvec: values at which bandwidth evaluated
% mse: vector of errors for each hvec
if nargin<9
    initvalues = zeros(3,1);
end

[nSubject,nVertex] = size(R);
if nSubject>nVertex 
    warning('N > V -- check that R is N x V')
end

ncrude = 10;
mindist = min(dist(dist~=0));
hvec = linspace(mindist-mindist/2,max(dist)/2,ncrude);
mse = zeros(ncrude,1);
K = 5;
famuniq = unique(familyID);
nfam = length(famuniq);
twinsmat = [MZtp1,MZtp2,DZtp1,DZtp2,MDti];

% use cross validation with FAMILY as the unit of observation
cvpart = cvpartition(nfam,'Kfold',K);
     
for j=1:ncrude
   
    h = hvec(j);
    index = (dist<=h);
    subR = R(:,index);
    wts = biweight(dist(index)/h)/h;
    sumwts = sum(wts);
    
    for k = 1:K
      % get indices for subjects
      fampart = famuniq(cvpart.training(k));
      subjIndex = ismember(familyID,fampart);
      tTwinsmat = twinsmat(subjIndex,:);
      tSubR = subR(subjIndex,:);
      n1 = sum(tTwinsmat(:,1)); n2 = sum(tTwinsmat(:,3)); n3 = sum(tTwinsmat(:,5));

      SM = tSubR(tTwinsmat(:,1),:).^2 + tSubR(tTwinsmat(:,2),:).^2; XM = 2*tSubR(tTwinsmat(:,1),:).*tSubR(tTwinsmat(:,2),:); % n1 x M
      SD = tSubR(tTwinsmat(:,3),:).^2 + tSubR(tTwinsmat(:,4),:).^2; XD = 2*tSubR(tTwinsmat(:,3),:).*tSubR(tTwinsmat(:,4),:); % n2 x M
      SMD = tSubR(tTwinsmat(:,5),:).^2; % n3 x M

      sumSM = sum(SM)*wts;
      sumXM = sum(XM)*wts;
      sumSD = sum(SD)*wts;
      sumXD = sum(XD)*wts;
      sumSMD = sum(SMD)*wts;

      options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton');
      [paramsA] = fminunc(@(params) wtdloglik_ACE(params, sumSM, sumXM, sumSD, sumXD, sumSMD, sumwts, n1, n2, n3), ...
           initvalues, options);     
      trainIndex =  logical(1 - subjIndex);
      tempR = R(trainIndex,dist==0);
      a = sum(exp(paramsA));
      mse(j) = mse(j) + mean((tempR.^2 - a).^2)/K;
    end
end

% fine search:

end

