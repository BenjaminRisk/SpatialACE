function [gradXa,gradXc] = gradient_psd(Xa,Xc,smSastar,smScstar,KJK)
%INPUT:
% Xa: V x Ka, where Ka is the rank of sigmaA
% Xc: V x Kc


% [~,~,unkernmat] = createkernmat(lat,long,h,false);
% smSastar = unkernmat*(2*S0+2*S1+S2)*unkernmat;
% KJK = unkernmat*ones(nVertex,nVertex)*unkernmat;

Sigmaa = Xa*Xa';
Sigmac = Xc*Xc';
  
     smA = smSastar - (4.5*Sigmaa + 5*Sigmac).*KJK;
     smC = smScstar - (5*Sigmaa + 6*Sigmac).*KJK;
 
     gradXa = -2*smA*Xa;
     gradXc = -2*smC*Xc;
   
end
