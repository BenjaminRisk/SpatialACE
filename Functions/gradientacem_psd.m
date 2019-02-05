function [gradXa,gradXc,gradXe] = gradientacem_psd(Xa,Xc,Xe,smSastar,smScstar,smS0,KJK)
%INPUT:
% Xa: V x Ka, where Ka is the rank of sigmaA
% Xc: V x Kc


% [~,~,unkernmat] = createkernmat(lat,long,h,false);
% smSastar = unkernmat*(2*S0+2*S1+S2)*unkernmat;
% KJK = unkernmat*ones(nVertex,nVertex)*unkernmat;


Sigmaa = Xa*Xa';
Sigmac = Xc*Xc';
Sigmae = Xe*Xe';  
     smA = smSastar - (4.5*Sigmaa + 5*Sigmac + 2*Sigmae).*KJK;
     smC = smScstar - (5*Sigmaa + 6*Sigmac + 2*Sigmae).*KJK;
     smE = 2*smS0 - 2*(Sigmaa + Sigmac + Sigmae).*KJK;
   
     gradXa = -2*smA*Xa;
     gradXc = -2*smC*Xc;
     gradXe = -2*smE*Xe;
end
