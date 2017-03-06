function [out] = fullcovacem_sandwich(R,sigmaem,familyst,lat,long,h,tol,outfull,bwonly)
% Benjamin Risk
% 6 March 2017
%FULLCOVACEM_SANDWICH calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive genetic, common environmental, and unique environmental variance components. 
%   This function is memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations and thus will not work for massive datasets. 
%
%INPUT:
% R: N x V matrix; residual matrix
% sigmaem: estimated MEASUREMENT ERROR; used in correction of S0
% familyst.MZtp1
% familyst.MZtp2
% familyst.DZtp1
% familyst.DZtp2
% familyst. MDti
% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar or vector of bandwidth(s); uses GCV to select output
%   covariances if h is a vector
% tol does not do anything if outfull=false, but is kept in this 
% argument position for backwards compatability with my old programs
% outfull true/false indicating whether to output the psd matrix
%   if false, no truncation occurs
% bwonly true/false indicating whether to only output the gcv results;
%   if true, the eigendecomposition of the smSA and smSC is NOT calculated

%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svalSEg
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% 
% smSA_symm
% smSC_symm
% smSEg_symm
% smSA_psd
% smSC_psd
% smSEg_psd

if nargin<9
    bwonly=false;
end
    
if nargin<8
    outfull=true;
end

if nargin<7
    tol=1e-03;
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

% get initial estimate of the unique environmental variance:
S0 = R'*R./N;
S0 = S0 - diag(sigmaem);

%% Calculate SA and SC matrix for cross-validation:

SA = 2*(S1 - S2);
SC = 2*S2 - S1;
SEg = S0 - S1;

clear S0;

nbw = length(h);
out.mseSA = zeros(nbw,1);
out.mseSC = out.mseSA;
out.mseSEg = out.mseSA;
for k=1:nbw

    kernmat = createkernmat(lat,long,h(k),true);
  
    % approximate GCV using the smoothed covariance matrix:
    tracekernmat = trace(kernmat);
    denomcovgcv = (1-tracekernmat^2/nVertex^2)^2;
    
    smSA_symm = kernmat*SA*kernmat';
    smSC_symm = kernmat*SC*kernmat';
    smSEg_symm = kernmat*SEg*kernmat';
    
    temp_sa = (SA - smSA_symm).^2;
    temp_sc = (SC - smSC_symm).^2;
    temp_seg = (SEg - smSEg_symm).^2;
    out.mseSA(k) = mean(temp_sa(:))/denomcovgcv;
    out.mseSC(k) = mean(temp_sc(:))/denomcovgcv;
    out.mseSEg(k) = mean(temp_seg(:))/denomcovgcv;
    
    if k==1 || out.mseSA(k) < min(out.mseSA(1:k-1))
        out.smSA_symm = smSA_symm;
    end
    
    if k==1 ||out.mseSC(k) < min(out.mseSC(1:k-1))
        out.smSC_symm = smSC_symm;
    end
    
    if k==1 || out.mseSEg(k) < min(out.mseSEg(1:k-1))
        out.smSEg_symm = smSEg_symm;
    end
end

out.hvec = h;
out.hvecmin = zeros(1,3);
%[~,b] = min(out.mseS1S2);
[~,b] = min(out.mseSA);
out.hvecmin(1) = h(b);
[~,b] = min(out.mseSC);
out.hvecmin(2) = h(b);
[~,b] = min(out.mseSEg);
out.hvecmin(3) = h(b);

out.smSA_symm = (out.smSA_symm+out.smSA_symm')./2;
out.smSC_symm = (out.smSC_symm+out.smSC_symm')./2;
out.smSEg_symm = (out.smSEg_symm+out.smSEg_symm')./2;

if ~bwonly
    [out.vecSA,out.valSA] = eig_descend(out.smSA_symm);
    [out.vecSC,out.valSC] = eig_descend(out.smSC_symm);
    [out.vecSEg,out.valSEg] = eig_descend(out.smSEg_symm);
    if outfull
        tindexA = find(out.valSA>tol);
        out.smSA_psd = out.vecSA(:,tindexA)*diag(out.valSA(tindexA))*out.vecSA(:,tindexA)';
        tindexC = find(out.valSC>tol);
        out.smSC_psd = out.vecSC(:,tindexC)*diag(out.valSC(tindexC))*out.vecSC(:,tindexC)';
        tindexEg = find(out.valSEg>tol);
        out.smSEg_psd = out.vecSEg(:,tindexEg)*diag(out.valSEg(tindexEg))*out.vecSEg(:,tindexEg)';
        out.h2_symm = diag(out.smSA_symm)./(diag(out.smSA_symm)+diag(out.smSC_symm) + diag(out.smSEg_symm));
        out.h2_psd = diag(out.smSA_psd)./(diag(out.smSA_psd)+diag(out.smSC_psd)+diag(out.smSEg_psd));
    end 
end

end
