function [svecSA,valSA,svecSC,valSC, svecSCnull, valSCnull] = fullcovace_original(R,MZtp1,MZtp2,DZtp1,DZtp2,MDti,lat,long,h,calcnull)
%FULLCOVACE calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function is memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations and thus will not work for massive datasets. 
%
%INPUT:
% R: N x V matrix; residual matrix
% MZtp1
% MZtp2
% DZtp1
% DZtp2
% MDti
% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar; bandwidth
% calcnull: Calculate SigmaC under null that SigmaA=0

%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
% svecSCnull: V x d matrix of eigenvectors under null hypothesis that SA=0
% svalSCnull: d x d matrix of associated eigenvalues
% mse: named cell giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for sigmaa, sigmac, and sigmac_null (if
%   calcnull==TRUE)
if nargin==9
    calcnull = 0;
end
if calcnull==0
    svecSCnull = [];
    valSCnull = [];
end

% create raw covariance estimates:
[n,nVertex] = size(R);
n1 = sum(MZtp1);
n2 = sum(DZtp1);

R11 = R(MZtp1,:);
R12 = R(MZtp2,:);
R1 = R11'*R12;
S1 = (R1 + R1')./2/n1;
clear R11 R12 R1;

R21 = R(DZtp1,:);
R22 = R(DZtp2,:);
R2 = R21'*R22; 
S2 = (R2 + R2')./2/n2;
clear R21 R22 R2;

S0 = R'*R./n;

%% Calculate SA matrix and its EVD:
SA = S0 + S1 - 2*S2;

%correct the diagonal:
for v=1:nVertex
     SA(v,v) = 2*S1(v,v) - 2*S2(v,v);
end
%SA = SA + diag(diag(S1)-diag(S0));

% truncate to positive eigenvalues:
[vecSA,valSA] = truncpartialevd(SA,n);
%[vecSA,valSA] = truncevd(SA);

clear SA;

%% Calculate SC matrix and its EVD:
SC = 2*S2 - 0.5*S0 - 0.5*S1;
%correct the diagonal:
for v=1:nVertex
      SC(v,v) = 2*S2(v,v)-S1(v,v);
end
SC = SC + 0.5*diag(diag(S0)-diag(S1));


% calculate SC matrix under null:
if calcnull
    % calculate SC under null hypothesis that SA=0
    SCnull = (S0 + S1 + S2)/3;
%     for v=1:nVertex
%         SCnull(v,v) = (S1(v,v)+S2(v,v))/2;
%     end
    SCnull = SCnull+diag(diag(S1)/6 + diag(S2)/6 - diag(S0)/3);
    clear S2 S0 S1
    [vecSCnull, valSCnull] = truncpartialevd(SCnull,n);
    %[vecSCnull, valSCnull] = truncevd(SCnull);
    clear SCnull;
end

clear S2 S0 S1;

[vecSC,valSC] = truncpartialevd(SC,n);
%[vecSC,valSC] = truncevd(SC);
clear SC;


%% Smooth the eigenvectors:
% create kernel matrix:
kernmat = createkernmat(lat,long,h,false);
svecSA = kernmat*vecSA;
svecSC = kernmat*vecSC;
if calcnull
    svecSCnull = kernmat*vecSCnull;
end
end