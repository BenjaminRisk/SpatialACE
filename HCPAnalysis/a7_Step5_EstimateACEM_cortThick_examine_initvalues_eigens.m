%----------------------
% Estimate ranks of covariance functions in Spatial ACEM
% Creates Web Supplement Figure: Screeplots. 
%------------------------------------

addpath ~/Dropbox/SpatialACDE/Programs/Functions
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load ~/Dropbox/SpatialACDE/Data/subjectData_cortThickForCovfun.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/SigmaemCortThickACEM.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_bandwidthselection_sw.mat

swresults{1}.hvecminL

ntwins = sum(familyst.MZtp1)+sum(familyst.DZtp1);
n1n2n3 = ntwins+sum(familyst.MDti);
nSubject = size(dataMatL,2);

sum(isnan(swresults{1}.outfull_swLR.valSA))
sum(isnan(swresults{1}.outfull_swLR.valSC))
sum(isnan(swresults{1}.outfull_swLR.valSEg))


valSA = swresults{1}.outfull_swLR.valSA;

valSC = swresults{1}.outfull_swLR.valSC;

valSEg = swresults{1}.outfull_swLR.valSEg;


valSA = swresults{1}.outfull_swLR.valSA(~isnan(swresults{1}.outfull_swLR.valSA));

valSC = swresults{1}.outfull_swLR.valSC(~isnan(swresults{1}.outfull_swLR.valSC));

valSA = valSA(1:900);

valSC = valSC(1:900);

valSEg = valSEg(1:1000);
%<-----------------


valSA(ntwins)
valSA(ntwins+1)

valSC(ntwins)
valSC(ntwins+1)


valSEg(nSubject-sum(familyst.MZtp1))
valSEg(nSubject-sum(familyst.MZtp1)+1)

valSEg(end-20:end)



h = subplot(2,3,1);

plot(1:length(valSA),valSA,'-.bo');
line([ntwins,ntwins],[-200,600],'Color','red') % number of twins

title('a) Eigenvalues of \Sigma_a')

subplot(2,3,4);

plot(1:length(valSA),cumsum(valSA),'--bo');
line([ntwins,ntwins],[0,20000],'Color','red') % number of twins

title('d) Sum of \Sigma_a')

subplot(2,3,2)

plot(1:length(valSC),valSC,'--bo');
line([ntwins,ntwins],[-200,600],'Color','red') % number of twins

title('b) Eigenvalues of \Sigma_c')

subplot(2,3,5)

plot(1:length(valSC),cumsum(valSC),'--bo');
line([ntwins,ntwins],[0,20000],'Color','red') % number of twins

title('e) Sum of \Sigma_c')

subplot(2,3,3)
plot(1:length(valSEg),valSEg,'--bo')
line([nSubject-sum(familyst.MZtp1),nSubject-sum(familyst.MZtp1)],[-200,600],'Color','green') % number of twins
title('c) Eigenvalues of \Sigma_{e,G}')

subplot(2,3,6)
plot(1:length(valSEg),cumsum(valSEg),'--bo')
line([nSubject-sum(familyst.MZtp1),nSubject-sum(familyst.MZtp1)],[0,20000],'Color','green') % number of twins
title('f) Sum of \Sigma_{e,G}')

saveas(h,'~/Dropbox/SpatialACDE/Documents/Figures/Figures_1200SubjectRelease/HCP_ACEM_Eigenvalues.jpg')

