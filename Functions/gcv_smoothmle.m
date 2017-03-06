function [out] = gcv_smoothmle(mleresults,familyst,lat,long,hvec,calcblups,useparfor)
% Benjamin Risk
% 6 March 2017
%GCV_SMOOTHMLE:
%INPUT:
% mleresults
% familyst
% lat
% long
% hvec
% calcblups: if true, calculate blups; defaults to false
% useparfor: if true, uses parfor to calculate the kernel matrix; defaults
%               to true
%OUTPUT:
% structure with fields as defined in mle_ACE

if nargin<6
    calcblups = false;
end

if nargin<7
    useparfor = true;
end

sigmasqa = mleresults.sigmasqA';
sigmasqc = mleresults.sigmasqC';
sigmasqe = mleresults.sigmasqE';
beta = mleresults.beta;
nbeta = size(beta,1);
MZtp1 = familyst.MZtp1;
MZtp2 = familyst.MZtp2;
DZtp1 = familyst.DZtp1;
DZtp2 = familyst.DZtp2;
MDti = familyst.MDti;

% Smoothed MLE:
nVertex = length(sigmasqa);
nk = length(hvec);
msevarcomp = zeros(3,nk);
msebeta = zeros(nbeta,nk);
out.smbeta = beta;
for k=1:nk
    %mykernmat = createkernmat(lat,long,hvec(k),false);
    if useparfor
            mykernmat = createkernmat_parfor(lat,long,hvec(k));
    else 
            mykernmat = createkernmat(lat,long,hvec(k));
    end
    tracekernmat = trace(mykernmat);
    denom = (1-tracekernmat/nVertex).^2;
    tempsmsigmasqa = mykernmat*sigmasqa;
    tempsmsigmasqc = mykernmat*sigmasqc;
    tempsmsigmasqe = mykernmat*sigmasqe;
    
    tempsmbeta = beta*mykernmat;
    
    msevarcomp(1,k) = mean((sigmasqa - tempsmsigmasqa).^2)./denom;
    msevarcomp(2,k) = mean((sigmasqc - tempsmsigmasqc).^2)./denom;
    msevarcomp(3,k) = mean((sigmasqe - tempsmsigmasqe).^2)./denom;
    msebeta(:,k) = mean((tempsmbeta - beta).^2,2)./denom;
    if k==1 || msevarcomp(1,k)<min(msevarcomp(1,1:k-1))
        % take transpose for compatibility with other functions
        out.smsigmasqa = tempsmsigmasqa';
    end
    if k==1 || msevarcomp(2,k)<min(msevarcomp(2,1:k-1))
        out.smsigmasqc = tempsmsigmasqc';
    end
    if k==1 || msevarcomp(3,k) < min(msevarcomp(3,1:k-1))
        out.smsigmasqe = tempsmsigmasqe';
    end
    for l=1:nbeta
        if k==1 || msebeta(l,k)<min(msebeta(l,1:k-1))
            out.smbeta(l,:) = tempsmbeta(l,:);
        end
    end
end

out.msevarcomp = msevarcomp;
out.hvec = hvec;
[~,b] = min(msevarcomp,[],2);
out.hvecmin = hvec(b);

out.msebeta = msebeta;
[~,b] = min(msebeta,[],2);
out.hvecbetamin = hvec(b);

ytemp = mleresults.fixefresid+mleresults.yfixed;
out.smyfixed = mleresults.xMat*out.smbeta;
out.smresid = ytemp - out.smyfixed;
out.h2 = (out.smsigmasqa./(out.smsigmasqa + out.smsigmasqc + out.smsigmasqe))';

% Calculate BLUPs:
if calcblups
N = length(MZtp1);
out.ai = zeros(N,nVertex);
out.aij = zeros(N,nVertex);
out.ci = zeros(N,nVertex);

    for v=1:nVertex

        sigmasqA = out.smsigmasqa(v);
        sigmasqC = out.smsigmasqc(v);
        sigmasqE = out.smsigmasqe(v);

        % MZs:
        n1 = sum(MZtp1);
        tindex = logical(MZtp1+MZtp2);
        tempSigma = kron(eye(n1),inv(sigmasqA*ones(2)+sigmasqC*ones(2)+sigmasqE*eye(2)));
        % Edited on 7 January 2016: Uses smoothed residuals
        interblups = tempSigma*out.smresid(tindex,v);
        out.ai(tindex,v) = kron(eye(n1),sigmasqA*ones(2))*interblups;
        out.ci(tindex,v) = kron(eye(n1),sigmasqC*ones(2))*interblups;

        % DZs:
        n2 = sum(DZtp1);
        tindex = logical(DZtp1+DZtp2);
        tempSigma = kron(eye(n2),inv(0.5*sigmasqA*eye(2)+0.5*sigmasqA*ones(2)+sigmasqC*ones(2)+sigmasqE*eye(2)));
        interblups = tempSigma*out.smresid(tindex,v);
        out.aij(tindex,v) = kron(eye(n2),0.5*sigmasqA*eye(2))*interblups;
        out.ai(tindex,v) = kron(eye(n2),0.5*sigmasqA*ones(2))*interblups;
        out.ci(tindex,v) = kron(eye(n2),sigmasqC*ones(2))*interblups;

        % Singletons:
        %n3 = sum(MDti);
        %tempSigma = eye(n3)./(out.sigmasqA+out.sigmasqC+out.sigmasqE);
        interblups = out.smresid(MDti,v)./(sigmasqA+sigmasqC+sigmasqE);
        out.ai(MDti,v) = sigmasqA*interblups;
        out.ci(MDti,v) = sigmasqC*interblups;
    end

    % Edited on 7 January 2016: Changed to smoothed fixed effects
    out.yhat = out.smyfixed+out.ai+out.aij+out.ci;
end

end
