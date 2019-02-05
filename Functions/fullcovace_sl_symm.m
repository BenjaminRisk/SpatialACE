function [out] = fullcovace_sl_symm(R,familyst,lat,long,h,sigmasqe)
%FULLCOVACE_SL_SYMM calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function is memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations and thus will not work for massive datasets. 
%   This function uses SL's estimator. 
%   "symm" denotes a symmetric, but not psd, matrix is output
%INPUT:
% R: N x V matrix; residual matrix
% familyst: structure with
%   MZtp1
%   MZtp2
%   DZtp1
%   DZtp2
%   MDti
% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar or vector of bandwidth(s); uses GCV to select output
%   covariances if h is a vector
% sigmasqe: optional input: estimates of unique environmental variance for
% use in h2 calculation
%OUTPUT:
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% 
% smSA_symm symmetric estimator, not psd
% smSA_psd truncated eigenvalues for psd
% smSC_symm
% smSC_psd

% NOTE: for 1,002 locations and a single BW, this implementation is 3,000 times faster
% than FSEM_cov with 1.5 seconds versus 4,209 seconds



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

%% Calculate SA matrix and its EVD:
SA = S0 + S1 - 2*S2;

% SL's version:
SA = SA - diag(diag(SA));


% [vecSA,valSA] = eig(SA);
% %re-order from largest to smallest:
% [valSA,tindex] = sort(diag(valSA),'descend');
% vecSA = vecSA(:,tindex);
% valSA = diag(valSA);

%% Calculate SC matrix and its EVD:
SC = 2*S2 - 0.5*S0 - 0.5*S1;

% SL's version:
%change this back
SC = SC - diag(diag(SC));

% [vecSC,valSC] = eig(SC);
% %re-order from largest to smallest:
% [valSC,tindex] = sort(diag(valSC),'descend');
% vecSC = vecSC(:,tindex);
% valSC = diag(valSC);

nbw = length(h);

for k=1:nbw

    % create kernel matrix:
    [~,rowsums,unkernmat] = createkernmat(lat,long,h(k),false);
    %svecSA = unkernmat*vecSA;
    %svecSC = unkernmat*vecSC;
    
    bigW = rowsums*rowsums' - unkernmat*unkernmat;
    
    if sum(isnan(bigW(:))) || sum(isinf(bigW(:)))
        error('bandwidth is too small - no neighbors within bandwidth')
    end
      
    % approximate GCV using the smoothed covariance matrix:
    % equivalent to the usual GCV when the points are equidistant on a
    % sphere, i.e., all entries of bigW are equal
    tracekernmat = sum(diag(unkernmat)./sqrt(diag(bigW)));
  
    denomcovgcv = (1-tracekernmat^2/nVertex^2)^2;
    smSA = unkernmat*SA*unkernmat./bigW;
    out.mseCov.smSA(k) = mean((SA(:)-smSA(:)).^2)/denomcovgcv;

    smSC = unkernmat*SC*unkernmat./bigW;
    
    out.mseCov.smSC(k) = mean((SC(:)-smSC(:)).^2)/denomcovgcv;

    if k==1 || out.mseCov.smSA(k)< min(out.mseCov.smSA(1:k-1))
        out.smSA_symm = smSA;
    end
    
    if k==1 || out.mseCov.smSC(k) < min(out.mseCov.smSC(1:k-1))
        out.smSC_symm = smSC;
    end
end

[out.svecSA,out.svalSA] = eig_descend(out.smSA_symm);
tindex = (out.svalSA>eps);
out.smSA_psd = out.svecSA(:,tindex)*diag(out.svalSA(tindex))*out.svecSA(:,tindex)';

[out.svecSC,out.svalSC] = eig_descend(out.smSC_symm);
tindex = out.svalSC>eps;
out.smSC_psd = out.svecSC(:,tindex)*diag(out.svalSC(tindex))*out.svecSC(:,tindex)';


out.hvec = h;
[~,b] = min(out.mseCov.smSA);
out.hvecmin(1) = h(b);
[~,b] = min(out.mseCov.smSC);
out.hvecmin(2) = h(b);

% option to calculate heritability:
if nargin>5
    out.h2_symm = diag(out.smSA_symm)./(diag(out.smSA_symm)+diag(out.smSC_symm)+sigmasqe');
    out.h2_psd = diag(out.smSA_psd)./(diag(out.smSA_psd)+diag(out.smSC_psd)+sigmasqe');
end
end
