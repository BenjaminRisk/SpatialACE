function [out] = estsigmasqem_cv(R,sigmasqA,sigmasqC,sigmasqE,lat,long,h,K)
%ESTSIGMAEM 
%       Estimates the covariance matrix for pooled additive genetic, 
%       common environmental, unique environmental.
%
%INPUT:
% R: N x V matrix; residual matrix
% sigmasqA
% sigmasqC
% sigmasqE: estimated MEASUREMENT ERROR PLUS UNIQUE ENVIRONMENTAL VARIANCES
% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar or vector of bandwidth(s); uses GCV to select output
%   covariances if h is a vector
% K number of partitions for cv


%OUTPUT:
% mse: named structure giving the mean squared error from the GCV approximation to
%   leave-one-out cross-validation for smSA, smSC, and smSCnull (if
%   calcnull==TRUE)
% hvec
% hvecmin
% sigmasqem

% create raw covariance estimates:
[N,nVertex] = size(R);
if N>nVertex 
    warning('N > V -- check that R is N x V')
end

cvpart = cvpartition(N,'Kfold',K);
nbw = length(h);
mse = zeros(K,nbw);

if nbw>1
    for m=1:nbw
      % create kernel matrix:
      [tempkern,rowsums,unkernmat] = createkernmat(lat,long,h(m),true);
      checkneighbors = any(diag(tempkern)>=0.999);
      if checkneighbors
          warning(['Bandwidth ' num2str(h(m)) ' is too small and ' num2str(checkneighbors) ' points have no neighbors'])
            mse(:,m) = NaN; 
      else

       for k=1:K
            trainR = R(cvpart.training(k),:);
            testR = (R(cvpart.test(k),:));
            nTrain = sum(cvpart.train(k));
            nTest = sum(cvpart.test(k));
            S0 = trainR'*trainR./nTrain;
            S0m = S0 - diag(diag(S0));
            clear S0;
            % huge communication overhead:
%             smSaceg = zeros(nVertexSub,1);
%             parfor v=1:nVertexSub
%                vunkernmat = unkernmat(:,subset(v));
%                bigW= rowsums(subset(v))^2 - vunkernmat'*vunkernmat;
%                smSaceg(v) = vunkernmat'*S0m*vunkernmat./bigW;
%             end

% many extra calculations but appears to be fastest!
            bigW = rowsums*rowsums' - unkernmat*unkernmat;
            smSaceg = (unkernmat*S0m*unkernmat)./bigW;
            
          %  allerror = (testRsq(:,subset) - repmat(smSaceg',nTest,1)).^2;
           % allerror = (testRsq - repmat(smSaceg',nTest,1)).^2;
            clear S0m bigW;
            %for j=1:nTest
            for j=1:1
                testS0 = testR(j,:)'*testR(j,:);
                allerror = (testS0 - smSaceg).^2;
                allerror = allerror - diag(diag(allerror)); % removes diagonal which is biased by measurement error
                mse(k,m) = mse(k,m) + sum(allerror(:))/(nVertex*(nVertex-1));
                clear testS0;
            end
       end
      end
    end


    out.mse = mean(mse);
    [~,b] = min(out.mse);
    out.hvecmin = h(b);
    out.hvec = h;
end

if nbw==1
    out.hvecmin=h;
end

S0 = R'*R./N;
S0 = S0 - diag(diag(S0));
[~,rowsums,unkernmat] = createkernmat(lat,long,out.hvecmin,true);
%smSaceg = zeros(nVertex,1); 
% parfor v=1:nVertex
%   vunkernmat = unkernmat(:,v);
%   bigW = rowsums(v)^2 - vunkernmat'*vunkernmat;
%   smSaceg(v) = vunkernmat'*S0*vunkernmat./bigW;
% end
bigW = diag(rowsums*rowsums' - unkernmat*unkernmat);
smSaceg = diag(unkernmat*S0*unkernmat);
clear unkernmat S0;
smSaceg = smSaceg./bigW;
inputSace = sigmasqA+sigmasqC+sigmasqE;
out.sigmasqem = inputSace - smSaceg';
end


