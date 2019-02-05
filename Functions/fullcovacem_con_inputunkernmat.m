function [ output] = fullcovacem_con_inputunkernmat(R,sigmasqem,dA,dC,dEg,familyst,unkernmat,maxiter,init,outfull)
%Constrained optimization of ACE covariance matrices
% R: N x V residual matrix
% sigmasqem: measurement error
% dA: rank of SigmaA
% dC: rank of SigmaC
% 


if nargin< 10
    outfull=1;
end

[N,nVertex] = size(R);

currentA = init.vecSA(:,1:dA)*diag(sqrt(init.valSA(1:dA)));
currentC = init.vecSC(:,1:dC)*diag(sqrt(init.valSC(1:dC)));
currentEg = init.vecSEg(:,1:dEg)*diag(sqrt(init.valSEg(1:dEg)));

% NOTE: here, the same bandwidth is used for all covariance functions.
% This could be revised in the future
KJK = (unkernmat*ones(nVertex,1))*(ones(1,nVertex)*unkernmat);

%--------------------------------->
% objects that do not change in gradient updates:
n1 = sum(familyst.MZtp1);
n2 = sum(familyst.DZtp1);

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
S0 = S0 - diag(sigmasqem);



Sastar = 2*S0+2*S1+S2;
clear S1;

Scstar = Sastar+S2;
clear S2;

smSastar = unkernmat*Sastar*unkernmat;
clear Sastar;
smScstar = unkernmat*Scstar*unkernmat;
clear Scstar;
% Updates for ACEM:
smS0 = unkernmat*S0*unkernmat;
clear unkernmat;

%<--------------------------------------   

% for simulations: lambda = 0.0001

% for hcp data:
%lambda = 0.00025;

lambdatol = 1e-9;
 lambda = 0.1;
 
 %lambda = 0.0005;% if too big, convergence problems
%lambda = 0.0001;

% algorithm has 4 V x V matrices:
% Sastar, Scstar, unkernmat, KJK
t=1;
propinitgrad = 1;
[gradA,gradC,gradEg] = gradientacem_psd(currentA,currentC,currentEg,smSastar,smScstar,smS0,KJK);
oldgradA = gradA;
oldgradC = gradC;
oldgradEg = gradEg;
newdeltaA = norm(gradA,'fro'); 
newdeltaC = norm(gradC,'fro');
newdeltaEg = norm(gradEg,'fro');
newnorm = sqrt(newdeltaA^2+newdeltaC^2+newdeltaEg^2);
initgradientnorm = newnorm;
%while propinitgrad>0.001 && t<=maxiter
counter=0;
while propinitgrad>0.0001 && t<=maxiter && lambda>lambdatol
    
    if t>1 && (newnorm >= gradnorm(t-1)) && counter<7
        lambda = lambda/2; 
        fprintf(['********************\n Norm of gradient increased. \n Decreasing lambda to ' num2str(lambda) '\n ********************\n' ]);
        currentA = oldA - lambda*oldgradA;
        currentC = oldC - lambda*oldgradC;
        currentEg = oldEg - lambda*oldgradEg;
%         [U,D] = svd(currentA);
%         currentA = U*D;
%         [U,D] = svd(currentC);
%     	currentC = U*D;
        [gradA,gradC,gradEg] = gradientacem_psd(currentA,currentC,currentEg,smSastar,smScstar,smS0,KJK);
        newdeltaA = norm(gradA,'fro'); 
        newdeltaC = norm(gradC,'fro');
        newdeltaEg = norm(gradEg,'fro');
        newnorm = sqrt(newdeltaA^2+newdeltaC^2+newdeltaEg^2); 
        counter=counter+1; %move after two failures
    else
       % t==1 || newnorm<gradnorm(t-1) || counter>=3
        counter=0;
        gradnorm(t) = newnorm;
        propinitgrad = gradnorm(t)/initgradientnorm;
        fprintf(['\nNorm of gradient: ' num2str(gradnorm(t)) '\n']);
        t=t+1;
        oldA = currentA;
        oldC = currentC;
        oldEg = currentEg;
        currentA = oldA - lambda*gradA;
        currentC = oldC - lambda*gradC;
        currentEg = oldEg - lambda*gradEg;
        oldgradA = gradA;
        oldgradC = gradC;
        oldgradEg = gradEg;
        [gradA,gradC,gradEg] = gradientacem_psd(currentA,currentC,currentEg,smSastar,smScstar,smS0,KJK);
        newdeltaA = norm(gradA,'fro'); 
        newdeltaC = norm(gradC,'fro');
        newdeltaEg = norm(gradEg,'fro');
        newnorm = sqrt(newdeltaA^2+newdeltaC^2+newdeltaEg^2);
    end
    
  
end

if t>=maxiter
    fprintf('***********************\n MaxIter reached \n');
    output.convergence = -1;

elseif lambda<=lambdatol
    fprintf(['*********\n Warning: did not converge. Current gradient norm: ' num2str(gradnorm(t-1)) '\n'...
        'Size relative to initial gradient: ' num2str(propinitgrad) '\n'...
        'If gradient is still large, try increasing dA and dC \n']);
    output.convergence = -2;
else
    output.convergence = 1;
end

     
output.gradnorm = gradnorm;
output.Xa = currentA;
output.Xc = currentC;
output.Xeg = currentEg;
if outfull
    output.smSA_psd = currentA*currentA';
    output.smSC_psd = currentC*currentC';
    output.smSEg_psd = currentEg*currentEg';
    output.h2 = diag(output.smSA_psd)./(diag(output.smSA_psd)+diag(output.smSC_psd)+diag(output.smSEg_psd));
end
end

