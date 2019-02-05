function [ output] = fullcovace_con_inputunkernmat(R,sigmasqe,dA,dC,familyst,unkernmat,maxiter,init,outfull)
%Constrained optimization of ACE covariance matrices
% R: N x V residual matrix
% dA: rank of SigmaA
% dC: rank of SigmaC
% familyst

% outfull: flag for whether to output the full covariance matrix

  
[nSubject,nVertex] = size(R);


% edit for big data:
if nargin< 9
     outfull=0;
end
% 
% [tvecA,tvalA] = eig_descend(init.smSA_psd);
% [tvecC,tvalC] = eig_descend(init.smSC_psd);
% 
% nzeigA = sum(diag(tvalA)>eps);
% nzeigC = sum(diag(tvalC)>eps);
% 
% 
% fprintf(['Note: inputed rank of SigmaA is ',num2str(dA),' \n and there are ',num2str(nzeigA) ' eigenvalues greater than eps\n'])
% fprintf(['\nNote: inputed rank of SigmaC is ',num2str(dC),' \n and there are ',num2str(nzeigC) ' eigenvalues greater than eps\n'])
% 
% newA = tvecA(:,1:dA)*sqrt(tvalA(1:dA,1:dA));
% newC = tvecC(:,1:dC)*sqrt(tvalC(1:dC,1:dC));
% 
% clear tvecA tvecC


KJK = unkernmat*ones(nVertex,nVertex)*unkernmat;

%--------------------------------->
% objects that do not change in gradient updates:
n1 = sum(familyst.MZtp1);
n2 = sum(familyst.DZtp1);
[N,nVertex] = size(R);

R11 = R(familyst.MZtp1,:);
R12 = R(familyst.MZtp2,:);
R1 = R11'*R12;
S1 = (R1 + R1')./2/n1;
clear R11 R12 R1;

S0 = R'*R./N;
S0 = S0 - diag(sigmasqe);
Sastar = 2*S0+2*S1;
clear S0 S1;

R21 = R(familyst.DZtp1,:);
R22 = R(familyst.DZtp2,:);
R2 = R21'*R22; 
S2 = (R2 + R2')./2/n2;
clear R21 R22 R2;

Sastar = Sastar+S2;

Scstar = Sastar+S2;
clear S2;

smSastar = unkernmat*Sastar*unkernmat;
clear Sastar;
smScstar = unkernmat*Scstar*unkernmat;
clear Scstar unkernmat; 
%<--------------------------------------   

% for simulations: lambda = 0.0001

% for hcp data:
%lambda = 0.00025;

% changed on 9 January 2017:
% Changed back to 0.001 on 10 January:
lambda = 0.001;
lambdatol = 1e-9;
%lambda = 0.01;

% algorithm has 4 V x V matrices:
% Sastar, Scstar, unkernmat, KJK
t=1;
propinitgrad = 1;
currentA = init.vecSA(:,1:dA)*diag(sqrt(init.valSA(1:dA)));
currentC = init.vecSC(:,1:dC)*diag(sqrt(init.valSC(1:dC)));

[gradA,gradC] = gradient_psd(currentA,currentC,smSastar,smScstar,KJK);
oldgradA = gradA;
oldgradC = gradC;
newdeltaA = norm(gradA,'fro'); 
newdeltaC = norm(gradC,'fro');
newnorm = sqrt(newdeltaA^2+newdeltaC^2);
initgradientnorm = newnorm;
counter=0;
while propinitgrad>0.0001 && t<=maxiter && lambda>lambdatol
    
    if t>1 && (newnorm >= gradnorm(t-1))
        lambda = lambda/2; 
        fprintf(['********************\n Norm of gradient increased. \n Decreasing lambda to ' num2str(lambda) '\n ********************\n' ]);
        currentA = oldA - lambda*oldgradA;
        currentC = oldC - lambda*oldgradC;
%         [U,D] = svd(currentA);
%         currentA = U*D;
%         [U,D] = svd(currentC);
%     	currentC = U*D;
        [gradA,gradC] = gradient_psd(currentA,currentC,smSastar,smScstar,KJK);
        newdeltaA = norm(gradA,'fro'); 
        newdeltaC = norm(gradC,'fro');
        newnorm = sqrt(newdeltaA^2+newdeltaC^2); 
        counter=counter+1; %move after four failures
    end
    
    if t==1 || newnorm<gradnorm(t-1) || counter>3
        counter=0;
        gradnorm(t) = newnorm;
        propinitgrad = gradnorm(t)/initgradientnorm;
        fprintf(['\nNorm of gradient: ' num2str(gradnorm(t)) '\n']);
        t=t+1;
        oldA = currentA;
        oldC = currentC;
        currentA = oldA - lambda*gradA;
        currentC = oldC - lambda*gradC;
        oldgradA = gradA;
        oldgradC = gradC;
        [gradA,gradC] = gradient_psd(currentA,currentC,smSastar,smScstar,KJK);
        newdeltaA = norm(gradA,'fro'); 
        newdeltaC = norm(gradC,'fro');
        newnorm = sqrt(newdeltaA^2+newdeltaC^2);
    end
end

if t>=maxiter
    fprintf('***********************\n MaxIter reached \n');
    output.convergence = -1;

elseif lambda<=lambdatol
    fprintf('*********\n Warning: did not converge. Initialize algorithm with smaller lambda');
    output.convergence = -2;
else
    output.converge = 1;
end



output.gradnorm = gradnorm;
output.Xa = currentA;
output.Xc = currentC;
if outfull
    output.smSA_psd = currentA*currentA';
    output.smSC_psd = currentC*currentC';
end
end

