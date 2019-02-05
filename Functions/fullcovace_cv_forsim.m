function [out] = fullcovace_cv_forsim(R,MZtp1,MZtp2,DZtp1,DZtp2,MDti,lat,long,h,sigmasqe,yfixed,calcnull)
%FULLCOVACE calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function is memory intensive because it constructs the full V x V residual matrix for intermediate 
%   calculations and thus will not work for massive datasets. 
%   Truncates eigenvalues BEFORE smoothing.
%INPUT:
% R: N x V matrix; residual matrix
% MZtp1
% MZtp2
% DZtp1
% DZtp2
% MDti
% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar or vector of bandwidth(s); uses GCV to select output
%   covariances if h is a vector
% calcnull: Calculate SigmaC under null that SigmaA=0

%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
% svecSCnull: V x d matrix of eigenvectors under null hypothesis that SA=0
% svalSCnull: d x d matrix of associated eigenvalues
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% 
% smSA
% smSC
% smSCnull

if nargin==9
    calcnull = 0;
end

% create raw covariance estimates:
[N,nVertex] = size(R);
n1 = sum(MZtp1);
n2 = sum(DZtp1);
n3 = sum(MDti);

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

S0 = R'*R./N;

%% Calculate SA matrix and its EVD:
SA = S0 + S1 - 2*S2;

SA = SA + diag(diag(S1)-diag(S0));

% truncate to positive eigenvalues:
%[vecSA,valSA] = truncpartialevd(SA,2*N);
[vecSA,valSA] = truncevd(SA);
psdSA = vecSA*valSA*vecSA';

%clear SA;

%% Calculate SC matrix and its EVD:
SC = 2*S2 - 0.5*S0 - 0.5*S1;
% correct the diagonal:
SC = SC + 0.5*diag(diag(S0)-diag(S1));


% calculate SC matrix under null:
if calcnull
    % calculate SC under null hypothesis that SA=0
    SCnull = (S0 + S1 + S2)/3;
    SCnull = SCnull+diag(diag(S1)/6 + diag(S2)/6 - diag(S0)/3);
    clear S2 S0 S1
    %[vecSCnull, valSCnull] = truncpartialevd(SCnull,2*N);
    [vecSCnull, valSCnull] = truncevd(SCnull);
    psdSCnull = vecSCnull*valSCnull*vecSCnull';
    %clear SCnull;
end

clear S2 S0 S1;

%[vecSC,valSC] = truncpartialevd(SC,2*N);
[vecSC,valSC] = truncevd(SC);
%clear SC;


psdSC = vecSC*valSC*vecSC';

nbw = length(h);
for k=1:nbw

    %% Smooth the eigenvectors:
    % create kernel matrix:
    kernmat = createkernmat(lat,long,h(k),false);
    svecSA = kernmat*vecSA;
    svecSC = kernmat*vecSC;
    
    % GCV using the smoothed covariance matrix:
    tracekernmat = trace(kernmat);
    denomcovgcv = (1-tracekernmat^2/nVertex^2)^2;
    smSA = svecSA*valSA*svecSA';
    out.mseCov.smSA(k) = mean((SA(:)-smSA(:)).^2)/denomcovgcv;

    smSC = svecSC*valSC*svecSC';
    out.mseCov.smSC(k) = mean((SC(:)-smSC(:)).^2)/denomcovgcv;

    if k==1 || out.mseCov.smSA(k)< min(out.mseCov.smSA(1:k-1))
        out.svecSA = svecSA;
        out.valSA = valSA;
        out.smSA = smSA;
    end
    
    if k==1 || out.mseCov.smSC(k) < min(out.mseCov.smSC(1:k-1))
        out.svecSC = svecSC;
        out.valSC = valSC;
        out.smSC = smSC;
    end
    
    if calcnull
        svecSCnull = kernmat*vecSCnull;
        smSCnull = svecSCnull*valSCnull*svecSCnull';
        out.mseCov.smSCnull(k) = mean((SCnull(:)-smSCnull(:)).^2)/denomcovgcv;  
        if k==1 || out.mseCov.smSCnull(k) < min(out.mseCov.smSCnull(1:k-1))
            out.valSCnull = valSCnull;
            out.svecSCnull = svecSCnull;
            out.smSCnull = smSCnull;
        end
    end
end

% Calculate eBLUPs:
out.ai = zeros(N,nVertex);
out.aij = zeros(N,nVertex);
out.ci = zeros(N,nVertex);
out.yhat = zeros(N,nVertex);

sigmaaJ2 = kron(smSA,ones(2,2));
sigmaaI2 = kron(smSA,eye(2));
sigmacJ2 = kron(smSC,ones(2,2));
sigmaeI2 = kron(diag(sigmasqe),eye(2));
invsigmaMZ = inv(sigmaaJ2 + sigmacJ2 + sigmaeI2);
invsigmaDZ = inv(0.5*sigmaaJ2+0.5*sigmaaI2+sigmacJ2+sigmaeI2);
invsigmaMD = inv(smSA+smSC+diag(sigmasqe));
out.invsigmaMZ = invsigmaMZ;
out.invsigmaDZ = invsigmaDZ;
out.invsigmaMD = invsigmaMD;


for i = 1:n1
    index = 2*i-1;
    yi = R(index:index+1,:);
    yi = yi(:);
    intermat = invsigmaMZ*yi;
    tempai = sigmaaJ2*intermat;
    tempai = reshape(tempai,2,nVertex);
    out.ai(index:index+1,:) = tempai;
    
    tempci = sigmacJ2*intermat;
    tempci = reshape(tempci,2,nVertex);
    out.ci(index:index+1,:) = tempci;
end    

for i = n1+1:n1+n2
    index = 2*i-1;
    yi = R(index:index+1,:);
    yi = yi(:);
    intermat = invsigmaDZ*yi;
    tempaij = 0.5*sigmaaI2*intermat;
    tempaij = reshape(tempaij,2,nVertex);
    out.aij(index:index+1,:) = tempaij;
    
    tempai = 0.5*sigmaaJ2*intermat;
    tempai = reshape(tempai,2,nVertex);
    out.ai(index:index+1,:) = tempai;
    
    tempci = sigmacJ2*intermat;
    tempci = reshape(tempci,2,nVertex);
    out.ci(index:index+1,:) = tempci;
end

for i=N-n3+1:N
    yi = R(i,:)';
    intermat = invsigmaMD*yi;
    out.ai(i,:) = smSA*intermat;
    out.ci(i,:) = smSC*intermat;
end
out.yhat = out.ai+out.ci+out.aij+yfixed;
out.hvec = h;
[~,b] = min(out.mseCov.smSA);
out.hvecmin(1) = h(b);
[~,b] = min(out.mseCov.smSC);
out.hvecmin(2) = h(b);
end
