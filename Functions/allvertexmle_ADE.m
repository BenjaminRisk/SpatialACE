function [out] = allvertexmle_ADE(ydata,xMat,familyst,useparfor)

% yData:    N subjects x V locations
% xMat:     N subjects x p
% familyst: 
    %MZtp1 
    %MZtp2 
    %DZtp1 
    %DZtp2 
    %MDti 
%useparfor either true or false; defaults to true which uses parfor
if nargin<4
    useparfor = true;
end

[N,nVertex] = size(ydata);
beta = zeros(size(xMat,2),nVertex);
pvalueBeta = beta;
sigmasqA = zeros(1,nVertex);
sigmasqD = zeros(1,nVertex);
sigmasqE = zeros(1,nVertex);
% mixefresid = zeros(N,nVertex);
fixefresid = zeros(N,nVertex);
%lrtA = zeros(1,nVertex);
%pvalueA = zeros(1,nVertex);
eFlag_a = zeros(1,nVertex);
%eFlag_n = zeros(1,nVertex);
n2loglik = zeros(1,nVertex);
% ai = zeros(N,nVertex);
% aij = zeros(N,nVertex);
% ci = zeros(N,nVertex);
yfixed = zeros(N,nVertex);
yhat = zeros(N,nVertex);

capmat = (xMat'*xMat)\xMat';
if useparfor
    parfor v=1:nVertex
    %for v=1:nVertex    
        estOLS = capmat*ydata(:,v); % for initializing coefficients
        results = mle_ADE(ydata(:,v), xMat, familyst, [0;0;0;estOLS]);
        beta(:,v) = results.betas;
        sigmasqA(v) = results.sigmasqA;
        sigmasqD(v) = results.sigmasqD;
        sigmasqE(v) = results.sigmasqE;
%         mixefresid(:,v) = results.mixefresid;
        fixefresid(:,v) = results.fixefresid;
      %  lrtA(v) = results.lrt;
     %   pvalueA(v) = results.pvalue;
        pvalueBeta(:,v) = results.betas_pvalues;

        eFlag_a(v) = results.exitflag_a;
 %       eFlag_n(v) = results.exitflag_n;
        n2loglik(v) = results.f;
%         ai(:,v) = results.ai;
%         aij(:,v) = results.aij;
%         ci(:,v) = results.ci;
        yfixed(:,v) = results.yfixed;
 %       yhat(:,v) = results.yhat;
    end
else
    for v=1:nVertex    
        estOLS = capmat*ydata(:,v); % for initializing coefficients
        results = mle_ADE(ydata(:,v), xMat, familyst, [0;0;0;estOLS]);
        beta(:,v) = results.betas;
        sigmasqA(v) = results.sigmasqA;
        sigmasqD(v) = results.sigmasqD;
        sigmasqE(v) = results.sigmasqE;
%         mixefresid(:,v) = results.mixefresid;
        fixefresid(:,v) = results.fixefresid;
%        lrtA(v) = results.lrt;
  %      pvalueA(v) = results.pvalue;
        pvalueBeta(:,v) = results.betas_pvalues;

        eFlag_a(v) = results.exitflag_a;
%        eFlag_n(v) = results.exitflag_n;
        n2loglik(v) = results.f;
%         ai(:,v) = results.ai;
%         aij(:,v) = results.aij;
%         ci(:,v) = results.ci;
        yfixed(:,v) = results.yfixed;
 %       yhat(:,v) = results.yhat;
    end
end

% annoying programming for parfor compatibility
out.beta = beta;
out.sigmasqA = sigmasqA;
out.sigmasqD = sigmasqD;
out.sigmasqE = sigmasqE;
%out.mixefresid = mixefresid;
out.fixefresid = fixefresid;
%out.lrtA = lrtA;
%out.pvalueA = pvalueA;
out.pvalueBeta = pvalueBeta;
out.eFlag_a = eFlag_a;
%out.eFlag_n = eFlag_n;
out.n2loglik = n2loglik;
% out.ai = ai;
% out.aij = aij;
% out.ci = ci;
out.yfixed = yfixed;
out.yhat = yhat;
out.xMat = xMat;
out.H2 = (sigmasqA./(sigmasqA+sigmasqD+sigmasqE))';

end
