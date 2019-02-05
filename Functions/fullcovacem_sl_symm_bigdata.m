function [out] = fullcovacem_sl_symm_bigdata(R,familyst,lat,long,h)
%FULLCOVACEM calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function does NOT calculate the low-rank decompositions or
%   truncated PSD estimators, which is computationally less expensive,
%   although memory intensive.
%   This function is also memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations. 
%   This function uses SL's estimator.
%   Update 6 February 2017: ACEM also estimates the covariance function for
%       correlated unique environmental effects.
%   "symm" denotes a symmetric, but not psd, matrix is output
%
%   NOTE: This function only calculates the covariances for ONE hemisphere.
%
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
%OUTPUT:
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% 
% smSA_symm symmetric estimator, not psd
% smSC_symm
% smSC_psd
% smSEg_symm unique environmental variance 
% smSEg_psd
% NOTE: for 1,002 locations and a single BW, this implementation is 3,000 times faster
% than FSEM_cov with 1.5 seconds versus 4,209 seconds

% NOTE: Currently, only variances and heritabilities are output, while
% this program calculates all correlations, which makes this implementation
% involve a lot of extraneous calculations. 

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
%SA = S0 + S1 - 2*S2;
% edited 6 February 2017:
SA = 2*(S1 - S2);

% SL's version:
SA = SA - diag(diag(SA));


% [vecSA,valSA] = eig(SA);
% %re-order from largest to smallest:
% [valSA,tindex] = sort(diag(valSA),'descend');
% vecSA = vecSA(:,tindex);
% valSA = diag(valSA);

%% Calculate SC matrix and its EVD:
%SC = 2*S2 - 0.5*S0 - 0.5*S1;
%edited 6 February 2017:
SC = 2*S2 - S1;

% SL's version:
SC = SC - diag(diag(SC));

% [vecSC,valSC] = eig(SC);
% %re-order from largest to smallest:
% [valSC,tindex] = sort(diag(valSC),'descend');
% vecSC = vecSC(:,tindex);
% valSC = diag(valSC);

% added 6 February 2017:
% Unique environmental genetic covariance
SEg = S0-S1;
SEg = SEg - diag(diag(SEg));


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

    smSEg = unkernmat*SEg*unkernmat./bigW;
    out.mseCov.smSEg(k) = mean((SEg(:)-smSEg(:)).^2)/denomcovgcv;
    
    
 % ONLY OUTPUT HERITABILITIES AND VARIANCES:
 
    if k==1 || out.mseCov.smSA(k)< min(out.mseCov.smSA(1:k-1))
        %out.smSA_symm = smSA;
        out.smSA_symm_variance = diag(smSA);
    end
    
    if k==1 || out.mseCov.smSC(k) < min(out.mseCov.smSC(1:k-1))
        %out.smSC_symm = smSC;     
        out.smSC_symm_variance = diag(smSC);
    end
    
    if k==1 || out.mseCov.smSEg(k) < min(out.mseCov.smSEg(1:k-1))
        %out.smSEg_symm = smSEg;
        out.smSEg_symm_variance = diag(smSEg);
    end
end
% 
% [out.svecSA,out.svalSA] = eig_descend(out.smSA_symm);
% tindex = (out.svalSA>eps);
% out.smSA_psd = out.svecSA(:,tindex)*diag(out.svalSA(tindex))*out.svecSA(:,tindex)';
% 
% [out.svecSC,out.svalSC] = eig_descend(out.smSC_symm);
% tindex = out.svalSC>eps;
% out.smSC_psd = out.svecSC(:,tindex)*diag(out.svalSC(tindex))*out.svecSC(:,tindex)';
% 
% [out.svecSEg,out.svalSEg] = eig_descend(out.smSEg_symm);
% tindex = out.svalSEg>eps;
% out.smSEg_psd = out.svecSEg(:,tindex)*diag(out.svalSEg(tindex))*out.svecSEg(:,tindex)';

out.hvec = h;
[~,b] = min(out.mseCov.smSA);
out.hvecmin(1) = h(b);
[~,b] = min(out.mseCov.smSC);
out.hvecmin(2) = h(b);
[~,b] = min(out.mseCov.smSEg);
out.hvecmin(3) = h(b);
out.h2_symm = out.smSA_symm_variance./(out.smSA_symm_variance+out.smSC_symm_variance+out.smSEg_symm_variance);
%out.h2_psd = diag(out.smSA_psd)./(diag(out.smSA_psd)+diag(out.smSC_psd)+diag(out.smSEg_psd));
end
