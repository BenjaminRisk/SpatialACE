function [out] = fullcovace_sandwich(R,sigmae,familyst,lat,long,h,tol,outfull,bwonly)
%FULLCOVACE_SANDWICH calculates a decomposition of the smoothed VxV covariance matrices
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
if nargin<9
    bwonly=false;
end
    
if nargin<8
    outfull=true;
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

nbw = length(h);
out.mseSASC = zeros(nbw,1);
for k=1:nbw

     kernmat = createkernmat(lat,long,h(k),true);
     
     
    % approximate GCV using the smoothed covariance matrix:
    tracekernmat = trace(kernmat);
    denomcovgcv = 2*(1-tracekernmat^2/nVertex^2)^2;
 
   
    smSA_symm = kernmat*SA*kernmat';
    smSC_symm = kernmat*SC*kernmat';
    
    % OLD APPROACH:
%     temp_s1 = (S1 - smSA_symm - smSC_symm).^2;
%     temp_s2 = (S2 - 0.5*smSA_symm - smSC_symm).^2;
%     out.mseS1S2(k) = (mean(temp_s1(:))+mean(temp_s2(:)))/denomcovgcv;
%     if k==1 || out.mseS1S2(k) < min(out.mseS1S2(1:k-1))
%         out.smSA_symm = smSA_symm;
%         out.smSC_symm = smSC_symm;
%     end
    
    % Edits 16 December 2016:
    temp_sa = (SA - smSA_symm).^2;
    temp_sc = (SC - smSC_symm).^2;
    out.mseSASC(k) = (mean(temp_sa(:))/2+mean(temp_sc(:))/2)/denomcovgcv;
    
    
    if k==1 || out.mseSASC(k) < min(out.mseSASC(1:k-1))
        out.smSA_symm = smSA_symm;
        out.smSC_symm = smSC_symm;
    end
end

out.hvec = h;
%[~,b] = min(out.mseS1S2);
[~,b] = min(out.mseSASC);
out.hvecmin = h(b);

if ~bwonly
    [out.vecSA,out.valSA] = eig_descend(out.smSA_symm);
    [out.vecSC,out.valSC] = eig_descend(out.smSC_symm);
    if outfull
        tindexA = find(out.valSA>tol);
        out.smSA_psd = out.vecSA(:,tindexA)*diag(out.valSA(tindexA))*out.vecSA(:,tindexA)';
        tindexC = find(out.valSC>tol);
        out.smSC_psd = out.vecSC(:,tindexC)*diag(out.valSC(tindexC))*out.vecSC(:,tindexC)';
        
        out.h2_symm = diag(out.smSA_symm)./(diag(out.smSA_symm)+diag(out.smSC_symm) + sigmae');
        out.h2_psd = diag(out.smSA_psd)./(diag(out.smSA_psd)+diag(out.smSC_psd)+sigmae');
    end 
end

end
