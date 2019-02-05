function [out] = fullcovace_sandwich_inputkernmat(R,sigmae,familyst,kernmat,tol,outfull,neigs)
%FULLCOVACE_SANDWICH_INPUTUNKERNMAT calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function is memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations and thus will not work for massive datasets. 
%   6 July 2016: Note: With the updated weighting approach, we can't smooth
%   post truncation. 
%INPUT:
% R: N x V matrix; residual matrix
% sigmae: estimated measurement error from mle; used in correction of S0
% familyst.MZtp1
% familyst.MZtp2
% familyst.DZtp1
% familyst.DZtp2
% familyst. MDti
% kernmat: V x V kernel matrix
% tol does not do anything if outfull=false, but is kept in this 
% argument position for backwards compatability with my old programs
% outfull true/false indicating whether to output the psd matrix
%   if false, no truncation occurs
% useeigs: integer indicating how many eigenvales to extract; the same
%   value is used for SA and SC
%

%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% 
% smSA_symm
% smSC_symm
% smSA_psd
% smSC_psd
if nargin<7
    neigs=0;
end
if nargin<6
    outfull=true;
end
if nargin<5
    tol = 1e-04;
end


% create raw covariance estimates:
[N,nVertex] = size(R);
n1 = sum(familyst.MZtp1);
n2 = sum(familyst.DZtp1);
n3 = sum(familyst.MDti);

R11 = R(familyst.MZtp1,:);
R12 = R(familyst.MZtp2,:);
R1 = R11'*R12;
S1 = (R1 + R1')./2/n1;
clear R11 R12 R1;

R21 = R(familyst.DZtp1,:);
R22 = R(familyst.DZtp2,:);
R2 = R21'*R22; 
S2 = (R2 + R2')./2/n2;
clear R21 R22 R2;

S0 = R'*R./N;
S0 = S0 - diag(sigmae);

%% Calculate SA and SC matrix for cross-validation:
%SA = S0plusS1 - 2*S2 +diag(diag(S1));
%SC = 2*S2 - 0.5*S0plusS1 - 0.5*diag(diag(S1));

SA = S0 + S1 - 2*S2;
SC = 2*S2 - 0.5*(S0+S1);

clear S0;

% approximate GCV using the smoothed covariance matrix:
%tracekernmat = trace(kernmat);
%denomcovgcv = (1-tracekernmat^2/nVertex^2)^2;
 
%out.mseCov.smSA(k) = mean((SA(:)-smSA(:)).^2)/denomcovgcv;
%out.mseCov.smSC(k) = mean((SC(:)-smSC(:)).^2)/denomcovgcv;
smSA_symm = kernmat*SA*kernmat';
smSC_symm = kernmat*SC*kernmat';
%temp_s1 = (S1 - out.smSA_symm - out.smSC_symm).^2;
%temp_s2 = (S2 - 0.5*out.smSA_symm - out.smSC_symm).^2;
    
%out.mseS1S2 = (mean(temp_s1(:))+mean(temp_s2(:)))/denomcovgcv;

smSA_symm = (smSA_symm + smSA_symm')./2;
smSC_symm = (smSC_symm + smSC_symm')./2;
[out.vecSA,out.valSA] = eig_descend(smSA_symm,neigs);
[out.vecSC,out.valSC] = eig_descend(smSC_symm,neigs);

if outfull
    tindexA = find(out.valSA>tol);
    out.smSA_psd = out.vecSA(:,tindexA)*diag(out.valSA(tindexA))*out.vecSA(:,tindexA)';
    tindexC = find(out.valSC>tol);
    out.smSC_psd = out.vecSC(:,tindexC)*diag(out.valSC(tindexC))*out.vecSC(:,tindexC)';
    out.smSA_symm = smSA_symm;
    clear smSA_symm;
    out.smSC_symm = smSC_symm;
end
end
