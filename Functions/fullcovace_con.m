function [ output] = fullcovace_con(R,sigmasqe,dA,dC,familyst,lat,long,h,maxiter,init,outfull)
%Constrained optimization of ACE covariance matrices
% R: N x V residual matrix
% dA: rank of SigmaA
% dC: rank of SigmaC
% 

if nargin < 9
% initialize from the truncated estimator:
    init = fullcovace_sandwich(R,sigmasqe,familyst,lat,long,h);
    %init = fullcovace_br_symm_v2(R,familyst,lat,long,h);
end

if nargin< 11
    outfull=1;
end

[N,nVertex] = size(R);
  

% commented out on 10 January 2017
% [tvecA,tvalA] = eig_descend(init.smSA_psd);
% [tvecC,tvalC] = eig_descend(init.smSC_psd);
% 
% nzeigA = sum(tvalA>eps);
% nzeigC = sum(tvalC>eps);
% 
% 
% fprintf(['Note: inputed rank of SigmaA is ',num2str(dA),' \n and there are ',num2str(nzeigA) ' eigenvalues greater than eps\n'])
% fprintf(['\nNote: inputed rank of SigmaC is ',num2str(dC),' \n and there are ',num2str(nzeigC) ' eigenvalues greater than eps\n'])
% 
% newA = tvecA(:,1:dA)*diag(sqrt(tvalA(1:dA)));
% newC = tvecC(:,1:dC)*diag(sqrt(tvalC(1:dC)));
% 
% clear tvecA tvecC
% 

currentA = init.vecSA(:,1:dA)*diag(sqrt(init.valSA(1:dA)));
currentC = init.vecSC(:,1:dC)*diag(sqrt(init.valSC(1:dC)));

[~,~,unkernmat] = createkernmat(lat,long,h,false);
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
S0 = S0 - diag(sigmasqe);

Sastar = 2*S0+2*S1+S2;
clear S0 S1;

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
% 10 January 2017 -- SIGNIFICANT MODIFICATIONS made to algorithm
%           to speed up convergence
%lambda = 0.001;
lambdatol = 1e-9;
 lambda = 0.1;
%lambda = 0.0005;% if too big, convergence problems
%lambda = 0.0001;

% algorithm has 4 V x V matrices:
% Sastar, Scstar, unkernmat, KJK
t=1;
propinitgrad = 1;
[gradA,gradC] = gradient_psd(currentA,currentC,smSastar,smScstar,KJK);
oldgradA = gradA;
oldgradC = gradC;
newdeltaA = norm(gradA,'fro'); 
newdeltaC = norm(gradC,'fro');
newnorm = sqrt(newdeltaA^2+newdeltaC^2);
initgradientnorm = newnorm;
%while propinitgrad>0.001 && t<=maxiter
counter=0;
while propinitgrad>0.001 && t<=maxiter && lambda>lambdatol
    
    if t>1 && (newnorm >= gradnorm(t-1)) && counter<2
        lambda = lambda/2; 
 %       fprintf(['********************\n Norm of gradient increased. \n Decreasing lambda to ' num2str(lambda) '\n ********************\n' ]);
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
        counter=counter+1; %move after two failures
    else
       % t==1 || newnorm<gradnorm(t-1) || counter>=3
        counter=0;
        gradnorm(t) = newnorm;
        propinitgrad = gradnorm(t)/initgradientnorm;
   %    fprintf(['\nNorm of gradient: ' num2str(gradnorm(t)) '\n']);
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
    output.convergence = 1;
end

     
output.gradnorm = gradnorm;
output.Xa = currentA;
output.Xc = currentC;
if outfull
    output.smSA_psd = currentA*currentA';
    output.smSC_psd = currentC*currentC';
    output.h2 = diag(output.smSA_psd)./(diag(output.smSA_psd)+diag(output.smSC_psd)+sigmasqe');
end
end

