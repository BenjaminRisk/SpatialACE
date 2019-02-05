function [out] = fullcovacem_sandwich_inputkernmat(R,sigmaem,familyst,kernmat,neigsSA,neigsSC,neigsSEg,outfull,tol)
%FULLCOVACEM_SANDWICH calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive genetic, common environmental, and unique environmental variance components. 
%   This function is memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations and thus will not work for massive datasets. 
%   6 July 2016: Note: With the updated weighting approach, we can't smooth
%   post truncation. 
%INPUT:
% R: N x V matrix; residual matrix
% sigmaem: estimated MEASUREMENT ERROR; used in correction of S0
% familyst.MZtp1
% familyst.MZtp2
% familyst.DZtp1
% familyst.DZtp2
% familyst. MDti
% kernmat
% neigsSA: equal to FALSE, indicating all eigenvalues are extracted, or an
%   integer denoting the number of eigenvalues to extract using eigs
% neigsSC:
% neigsSEg:
%
%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svalSEg
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
%
% smSA_symm
% smSC_symm
% smSEg_symm
% smSA_psd
% smSC_psd
% smSEg_psd

if nargin<9
    tol=1e-04;
end
    
if nargin<8
    outfull=true;
end

if nargin<5
    neigsSA=0; %neigs = 0 is a switch that uses full EVD instead of eigs with number
                % of eigenvalues equals to neigs
    neigsSC = 0;
    neigsSEg = 0;
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
S0 = S0 - diag(sigmaem);

SA = 2*(S1 - S2);
SC = 2*S2 - S1;
clear S2;
SEg = S0 - S1;

clear S0 S1;

smSA_symm = kernmat*SA*kernmat';
clear SA;
smSC_symm = kernmat*SC*kernmat';
clear SC;
smSEg_symm = kernmat*SEg*kernmat';
clear SEg;    

smSA_symm = (smSA_symm + smSA_symm')./2;
smSC_symm = (smSC_symm + smSC_symm')./2;
smSEg_symm = (smSEg_symm + smSEg_symm')./2;

[out.vecSA,out.valSA] = eig_descend(smSA_symm,neigsSA);
[out.vecSC,out.valSC] = eig_descend(smSC_symm,neigsSC);
[out.vecSEg,out.valSEg] = eig_descend(smSEg_symm,neigsSEg);

if outfull
    out.smSA_symm = smSA_symm;
    out.smSC_symm = smSC_symm;
    out.smSEg_symm = smSEg_symm;
    tindexA = find(out.valSA>tol);
    out.smSA_psd = out.vecSA(:,tindexA)*diag(out.valSA(tindexA))*out.vecSA(:,tindexA)';
    tindexC = find(out.valSC>tol);
    out.smSC_psd = out.vecSC(:,tindexC)*diag(out.valSC(tindexC))*out.vecSC(:,tindexC)';
    tindexEg = find(out.valSEg>tol);
    out.smSEg_psd = out.vecSEg(:,tindexEg)*diag(out.valSEg(tindexEg))*out.vecSEg(:,tindexEg)';
end

end
