function [out] = estsigmasqem_gcv(R,familyst,sigmasqA,sigmasqC,sigmasqE,lat,long,h)
%ESTSIGMAEM 
%       Estimates the covariance matrix for pooled additive genetic, 
%       common environmental, unique environmental.
%
%INPUT:
% R: N x V matrix; residual matrix
% sigmasqA
% sigmasqC
% sigmasqE: estimated MEASUREMENT ERROR PLUS UNIQUE ENVIRONMENTAL VARIANCES
% lat
% long
% h

% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar or vector of bandwidth(s); uses GCV to select output
%   covariances if h is a vector

%OUTPUT:
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% hvec
% hvecmin
% sigmasqem

[N,nVertex] = size(R);
if N>nVertex 
    warning('N > V -- check that R is N x V')
end

n1 = sum(familyst.MZtp1);
nbw = length(h);
mse = zeros(1,nbw);
R11 = R(familyst.MZtp1,:);
R12 = R(familyst.MZtp2,:);
S1 = (R11'*R12 + R12'*R11)/2/n1;

if nbw>1
    for m=1:nbw
      % create kernel matrix:
        kernmat = createkernmat(lat,long,h(m),true);
        checkneighbors = any(diag(kernmat)>=0.999);
        if checkneighbors
              warning(['For bandwidth ' num2str(h(m)), ' ' num2str(checkneighbors) ' points have no neighbors'])
        end
                 % approximate GCV using the smoothed covariance matrix:
        tracekernmat = trace(kernmat);
        denomcovgcv = (1-tracekernmat^2/nVertex^2)^2;
        smSaceg = kernmat*S1*kernmat';
        temp_sa = norm(S1 - smSaceg,'fro')^2;
        mse(m) = temp_sa/nVertex^2/denomcovgcv;
    end
out.mse = mse;
[~,b] = min(mse);
out.hvecmin = h(b);
out.hvec = h;
end
clear S1 smSaceg;

if nbw==1
    out.hvecmin=h;
end

S0 = R'*R./N;
S0 = S0 - diag(diag(S0));
[~,rowsums,unkernmat] = createkernmat(lat,long,out.hvecmin,true);
% too much communication overhead for parfor
%smSaceg = zeros(nVertex,1); 
% parfor v=1:nVertex
%   vunkernmat = unkernmat(:,v);
%   bigW = rowsums(v)^2 - vunkernmat'*vunkernmat;
%   smSaceg(v) = vunkernmat'*S0*vunkernmat./bigW;
% end
bigW = rowsums*rowsums' - unkernmat*unkernmat;
smSaceg = diag(unkernmat*S0*unkernmat)./diag(bigW);
inputSace = sigmasqA+sigmasqC+sigmasqE;
out.sigmasqem = inputSace - smSaceg';
end


