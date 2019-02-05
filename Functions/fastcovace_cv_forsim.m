function [out] = fastcovace_cv_forsim(R,MZtp1,MZtp2,DZtp1,DZtp2,MDti,lat,long,hvec,sigmasqe,yfixed,lomem)
%FASTCOVACE calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function uses random projections to avoid calculating the V x V
%   residual matrices and thus is scalable to large covariance matrices.
%svecSA,valSA,svecSC,valSC,

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
% calcnull: calculate SigmaC under the null hypothesis that SigmaA = 0
% lomem: use sparse matrices; avoid explicitly calculating covariance matrix
%        slower but scalable to massive matrices

%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
% svecSCnull: V x d matrix of eigenvectors under null hypothesis that SA=0
% svalSCnull: d x d matrix of associated eigenvalues

[N,nVertex] = size(R);

% rank of random projection matrix:
%dimrp = N+20;
dimrp = min(3*N,nVertex/2);

% generate a matrix that is approximately orthogonal:
Q = randn(nVertex,dimrp)/sqrt(dimrp);
% OPTIONAL: orthogonalize
[Q, ~, ~] = svd(Q,'econ');

% reducing bias slightly increases the MSE:
Q = Q*sqrt(nVertex/dimrp);

n1 = sum(MZtp1);
n2 = sum(DZtp1);
n3 = sum(MDti);

projR = R*Q;
R11 = R(MZtp1,:);
R12 = R(MZtp2,:);
projR11 = projR(MZtp1,:);
projR12 = projR(MZtp2,:);
projR1 = projR11'*projR12;
projS1 = (projR1 + projR1')./2/n1;

diagS0 = sum(R.^2)/N;
diagS1 = sum(R11.*R12)/n1;
clear R11 R12 R1;

R21 = R(DZtp1,:);
R22 = R(DZtp2,:);
%diagS2 = sum(R21.*R22)/n2;

projR21 = projR(DZtp1,:);
projR22 = projR(DZtp2,:);
projR2 = projR21'*projR22;
projS2 = (projR2 + projR2')./2/n2;
clear R21 R22 R2;

projS0 = projR'*projR./N;

%% Calculate SA matrix and its EVD:
projSA = projS0 + projS1 - 2*projS2;

%correct the diagonal:
diagSA = sparse(1:nVertex,1:nVertex,diagS1 - diagS0);
projdiagSA = Q'*diagSA*Q;
projSA = projSA + projdiagSA;
% truncate to positive eigenvalues:
projSA = (projSA+projSA')/2;
[vecprojSA,valprojSA] = eig(projSA);
valprojSA = diag(valprojSA);
[~,sindex] = sort(valprojSA,'descend');
vecprojSA = vecprojSA(:,sindex);
valprojSA = valprojSA(sindex);
index = valprojSA > (0.01*valprojSA(1));
valSA = diag(valprojSA(index));
psdprojSA = vecprojSA(:,index)*valSA*vecprojSA(:,index)';
vecSAall = Q*vecprojSA;
vecSA = vecSAall(:,index);

%valSAall = diag(valprojSA);
if lomem==false
    SA = vecSAall*diag(valprojSA)*vecSAall';
    psdSA = vecSA*valSA*vecSA';
end

%% Calculate SC matrix and its EVD:
projSC = 2*projS2 - 0.5*projS0 - 0.5*projS1;

% correct the diagonal:
diagSC = 0.5*(diagS0 - diagS1);
diagSC = sparse(1:nVertex,1:nVertex,diagSC);
projdiagSC = Q'*diagSC*Q;
projSC = projSC+projdiagSC;
projSC = (projSC+projSC')/2;
[vecprojSC,valprojSC] = eig (projSC);
valprojSC = diag(valprojSC);
[~,sindex] = sort(valprojSC,'descend');
vecprojSC = vecprojSC(:,sindex);
valprojSC = valprojSC(sindex);
index = valprojSC> (0.01*valprojSC(1));
valSC = diag(valprojSC(index));
psdprojSC = vecprojSC(:,index)*valSC*vecprojSC(:,index)';
vecSCall = Q*vecprojSC;
vecSC = vecSCall(:,index);

%valSCall = diag(valprojSC);
if lomem==false
    SC = vecSCall*diag(valprojSC)*vecSCall';
    psdSC = vecSC*valSC*vecSC';
end


%NOTE: Consider treating KQQ' as the smoothing matrix in gcv
out.mseCov = zeros(1,length(hvec));
for k=1:length(hvec)
    kernmat = createkernmat(lat,long,hvec(k),lomem);
    svecSA = kernmat*vecSA;
    svecSC = kernmat*vecSC;
    %calculate elements of covariance matrix in blocks:
    % GCV using the smoothed covariance matrix:
    tracekernmat = trace(kernmat);
    % 1 - 2*tracekernmat^2/2/V^2
    denomcovgcv = (1-tracekernmat^2/nVertex^2)^2;
    
    if lomem==false
        smSA = svecSA*valSA*svecSA';
        smSC = svecSC*valSC*svecSC';
        out.mseCov(k) = mean(([SA(:);SC(:)]-[smSA(:);smSC(:)]).^2)/denomcovgcv;
    end
    
    if k==1 || out.mseCov(k)< min(out.mseCov(1:k-1))
        out.svecSA = svecSA;
        out.svecSC = svecSC;
        KQ = kernmat*Q;
        if lomem==false
            out.smSA = smSA;
            out.smSC = smSC;
        end
    end
    clear kernmat
end
clear smSA smSC
out.hvec = hvec;
[~,b] = min(out.mseCov);
out.hvecmin = hvec(b);

% Calculate eBLUPs:
out.ai = zeros(N,nVertex);
out.aij = zeros(N,nVertex);
out.ci = zeros(N,nVertex);

% Create low rank decompositions for prediction
projsigmaMD = psdprojSA+psdprojSC;
projsigmaDZ = 0.5*kron(psdprojSA,ones(2))+0.5*kron(psdprojSA,eye(2))+kron(psdprojSC,ones(2));
spsigmae = sparse(1:nVertex,1:nVertex,sigmasqe);


% for now, construct full inverse matrices
invsigmaMZ =  smw(kron(KQ,ones(2,1)),projsigmaMD,kron(spsigmae,eye(2)));
invsigmaDZ = smw(kron(KQ,eye(2)),projsigmaDZ,kron(spsigmae,eye(2))); %kron(KQ,speye(2)) much slower
invsigmaMD = smw(KQ,projsigmaMD,spsigmae);
out.invsigmaMZ = invsigmaMZ;
out.invsigmaDZ = invsigmaDZ;
out.invsigmaMD = invsigmaMD;

for i = 1:n1
    index = 2*i-1;
    yi = R(index:index+1,:);
    yi = yi(:);
    intermat = invsigmaMZ*yi;
    tempai = kron(out.smSA,ones(2))*intermat;
    tempai = reshape(tempai,2,nVertex);
    out.ai(index:index+1,:) = tempai;
    
    tempci = kron(out.smSC,ones(2))*intermat;
    tempci = reshape(tempci,2,nVertex);
    out.ci(index:index+1,:) = tempci;
end    

for i = n1+1:n1+n2
    index = 2*i-1;
    yi = R(index:index+1,:);
    yi = yi(:);
    intermat = invsigmaDZ*yi;
    tempaij = 0.5*kron(out.smSA,speye(2))*intermat;
    tempaij = reshape(tempaij,2,nVertex);
    out.aij(index:index+1,:) = tempaij;
    
    tempai = 0.5*kron(out.smSA,ones(2))*intermat;
    tempai = reshape(tempai,2,nVertex);
    out.ai(index:index+1,:) = tempai;
    
    tempci = kron(out.smSC,ones(2))*intermat;
    tempci = reshape(tempci,2,nVertex);
    out.ci(index:index+1,:) = tempci;
end

for i=N-n3+1:N
    yi = R(i,:)';
    intermat = invsigmaMD*yi;
    out.ai(i,:) = out.smSA*intermat;
    out.ci(i,:) = out.smSC*intermat;
end
out.yhat = out.ai+out.ci+out.aij+yfixed;
end