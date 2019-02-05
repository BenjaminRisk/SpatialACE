%-----------------------------------------
% Calculate heritability and variances from low-rank decompositions.
% Create figures. 
% Create CIFTI files. 
%-----------------------------------------

addpath ~/Dropbox/SpatialACDE/JASA_Reproducibility/Functions/

% cifti-matlab available at https://github.com/Washington-University/cifti-matlab
addpath ~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8


% the supporting files are available on my github:
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R

load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/mle_cortThick676_HIGH.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/mwle_CortThick676_HIGH.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat

load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/SigmaemCortThickACEM.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_bandwidthselection_sw.mat



nVertexL = length(latL);
nVertexR = length(latR);
nVertex = nVertexL + nVertexR;


%% create variance estimates from PSD-ACE:
% compare to fine tune:
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_finetune_K=1.mat
psdresults2 = psdresults;

psdresults2{1}.smsigmasqa = zeros(nVertex,1);
psdresults2{1}.smsigmasqc = zeros(nVertex,1);
psdresults2{1}.smsigmasqeg = zeros(nVertex,1);

for v=1:nVertex
    psdresults2{1}.smsigmasqa(v) = psdresults2{1}.outfull_psdLR.Xa(v,:)*psdresults2{1}.outfull_psdLR.Xa(v,:)';
    psdresults2{1}.smsigmasqc(v) = psdresults2{1}.outfull_psdLR.Xc(v,:)*psdresults2{1}.outfull_psdLR.Xc(v,:)';
    psdresults2{1}.smsigmasqeg(v) = psdresults2{1}.outfull_psdLR.Xeg(v,:)*psdresults2{1}.outfull_psdLR.Xeg(v,:)';
end
psdresults2{1}.h2 = psdresults2{1}.smsigmasqa./(psdresults2{1}.smsigmasqa + psdresults2{1}.smsigmasqc + psdresults2{1}.smsigmasqeg);
mean(psdresults2{1}.h2)
std(psdresults2{1}.h2)
max(psdresults2{1}.h2)
min(psdresults2{1}.h2)
input = [psdresults2{1}.smsigmasqa,psdresults2{1}.smsigmasqc,psdresults2{1}.smsigmasqeg,sigmaemresults{1}.smsigmasqemLR',...
    psdresults2{1}.h2];
outfile = '~/Dropbox/SpatialACDE/Results/heritability_PSDACEM_finetune'; 
myft_write_cifti(input,outfile,'~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii',...
    '~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L',...
    '~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R');



mean(input)

%% The results from the original estimate are used in figures and manuscript summaries:
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_divideandconquer_K=1.mat

psdresults{1}.smsigmasqa = zeros(nVertex,1);
psdresults{1}.smsigmasqc = zeros(nVertex,1);
psdresults{1}.smsigmasqeg = zeros(nVertex,1);

for v=1:nVertex
    psdresults{1}.smsigmasqa(v) = psdresults{1}.outfull_psdLR.Xa(v,:)*psdresults{1}.outfull_psdLR.Xa(v,:)';
    psdresults{1}.smsigmasqc(v) = psdresults{1}.outfull_psdLR.Xc(v,:)*psdresults{1}.outfull_psdLR.Xc(v,:)';
    psdresults{1}.smsigmasqeg(v) = psdresults{1}.outfull_psdLR.Xeg(v,:)*psdresults{1}.outfull_psdLR.Xeg(v,:)';
end
psdresults{1}.h2 = psdresults{1}.smsigmasqa./(psdresults{1}.smsigmasqa + psdresults{1}.smsigmasqc + psdresults{1}.smsigmasqeg);


input = [psdresults{1}.smsigmasqa,psdresults{1}.smsigmasqc,psdresults{1}.smsigmasqeg,sigmaemresults{1}.smsigmasqemLR',...
    psdresults{1}.h2];
outfile = '~/Dropbox/SpatialACDE/Results/heritability_PSDACEM'; 
myft_write_cifti(input,outfile);


mean(psdresults{1}.smsigmasqa)
mean(psdresults{1}.smsigmasqc)
mean(psdresults{1}.smsigmasqeg)
mean(sigmaemresults{1}.smsigmasqemLR)

mean(psdresults{1}.h2)
std(psdresults{1}.h2)
max(psdresults{1}.h2)
min(psdresults{1}.h2)

mean(abs(psdresults2{1}.h2 - psdresults{1}.h2))
max(abs(psdresults2{1}.h2 - psdresults{1}.h2))
sum(abs(psdresults2{1}.h2 - psdresults{1}.h2))
sum(psdresults2{1}.h2 > psdresults{1}.h2)
% the max change in heritability equaled 0.0048. We take this as
% evidence of adequate convergence. Note in the finetune data,
% the heritability actually increased in all locations, although
% negligible.

% other notes from .out files:
% norm of initial gradient: 151,260.
% norm of grad at end of
% AllEstimatesCortThickACEM_divideandconquer_K=1.mat: 29.744, or 0.02%
% AllEstimatesCortThickACEM_finetune_K=1: 15.22, or 0.01%


%% ASIDE: Tally number of times mle ACE was higher than mle ADE
mean([mleresultsADE_L_HIGH.n2loglik,mleresultsADE_R_HIGH.n2loglik]<[mleresultsL_HIGH.n2loglik,mleresultsR_HIGH.n2loglik])
sum([mleresultsADE_L_HIGH.n2loglik,mleresultsADE_R_HIGH.n2loglik]<[mleresultsL_HIGH.n2loglik,mleresultsR_HIGH.n2loglik])
length([mleresultsADE_L_HIGH.n2loglik,mleresultsADE_R_HIGH.n2loglik])

% create measurement error cifti
outfile = '~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/measurementError';
myft_write_cifti(sigmaemresults{1}.smsigmasqemLR',outfile);

% cross-validation results for smoothing mle estimates of variance
% components:
subplot(2,2,1)
plot(smmleL_HIGH.hvec,smmleL_HIGH.msevarcomp)
title('MLE \sigma_a^2, \sigma_c^2, \sigma_e^2')
legend({'\sigma_a^2','\sigma_c^2','\sigma_e^2'})


% cross-validation results for smoothing mle estimates of covariates:
subplot(2,2,2)
plot(smmleL_HIGH.hvec,smmleL_HIGH.msebeta)
ylim([0,6e-4])
title('MLE covariates')
smmleL_HIGH.hvecbetamin 
smmleR_HIGH.hvecbetamin

legend({'Intercept','Age','Gender','TIV'})



% generalized cross validation for additive genetic + common environmental + unique
% environmental
subplot(2,2,3)
plot(sigmaemresults{1}.estsmsigmasqemL.hvec,sigmaemresults{1}.estsmsigmasqemL.mse)
title('\Sigma_G')

subplot(2,2,4)
temp = [swresults{1}.mseSA_L,swresults{1}.mseSC_L,swresults{1}.mseSEg_L];
subplot(2,2,4)
plot([1.18,1.25,1.3,1.35,1.4,1.5,1.75],temp)
title('\Sigma_a   \Sigma_c   and \Sigma_{e,G}')
legend({'\Sigma_a','\Sigma_c','\Sigma_{e,G}'})

saveas(gca,'~/Dropbox/SpatialACDE/Documents/Figures/Figures_1200SubjectRelease/GCV_HCP.jpg')
swresults{1}.hvecminL


% look at weights for some of the chosen bandwidths:
nVertexL = length(latL);
indexx = randsample(nVertexL,5000);
kern1pt4 = createkernmat(latL,longL,1.4,true);
kernmatsort =  zeros(5000,nVertexL);

for v=1:5000
    kernmatsort(v,:) = sort(kern1pt4(indexx(v),:),'descend');
end
temp = mean(kernmatsort,1);
temp(1:10)


kern1pt3 = createkernmat(latL,longL,1.3,true);
kernmatsort1pt3 = zeros(5000,nVertexL);
for v=1:5000
    kernmatsort1pt3(v,:) = sort(kern1pt3(indexx(v),:),'descend');
end
temp = mean(kernmatsort1pt3,1);
temp(1:10)

%% Calculate heritability in SMLE, MLE, and MWLE
% Create CIFTI files with left and right hemis for smmle, mwle, mle, and conpsd:
smsigmasqa_smmle = [smmleL_HIGH.smsigmasqa,smmleR_HIGH.smsigmasqa];
smsigmasqc_smmle = [smmleL_HIGH.smsigmasqc,smmleR_HIGH.smsigmasqc];
smsigmasqe_smmle = [smmleL_HIGH.smsigmasqe,smmleR_HIGH.smsigmasqe];
h2_smmle = [smmleL_HIGH.h2',smmleR_HIGH.h2'];

mean(h2_smmle)
std(h2_smmle)

% to gain insight into higher heritability in PSD-ACE, correct
% h2_smmle with measurement error correction:
sigmag = smsigmasqa_smmle+smsigmasqc_smmle+smsigmasqe_smmle - sigmaemresults{1}.smsigmasqemLR;
sum(sigmag<0)
% none are negative:)
mean(sigmaemresults{1}.smsigmasqemLR)
h2_smmle_mec = smsigmasqa_smmle./sigmag;
mean(h2_smmle_mec)
std(h2_smmle_mec)


input = [smsigmasqa_smmle',smsigmasqc_smmle',smsigmasqe_smmle',h2_smmle'];
outfile = '~/Dropbox/SpatialACDE/Results/heritability_smmle';
myft_write_cifti(input,outfile);

sigmasqa_mle = [mleresultsL_HIGH.sigmasqA,mleresultsR_HIGH.sigmasqA];
sigmasqc_mle = [mleresultsL_HIGH.sigmasqC,mleresultsR_HIGH.sigmasqC];
sigmasqe_mle = [mleresultsL_HIGH.sigmasqE,mleresultsR_HIGH.sigmasqE];
h2_mle = [mleresultsL_HIGH.h2',mleresultsR_HIGH.h2'];
input = [sigmasqa_mle',sigmasqc_mle',sigmasqe_mle',h2_mle'];
outfile = '~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/heritability_mle';
myft_write_cifti(input,outfile);
mean(h2_mle)
std(h2_mle)
max(h2_mle)

mean(psdresults{1}.h2'-h2_mle)

sigmag_mle = sigmasqa_mle  + sigmasqc_mle + sigmasqe_mle -   sigmaemresults{1}.smsigmasqemLR;
h2_mle_mec = sigmasqa_mle./sigmag_mle;
mean(h2_mle_mec)
std(h2_mle_mec)


subplot(2,3,1)
hist(psdresults{1}.smsigmasqa,20)
subplot(2,3,2)
hist(psdresults{1}.smsigmasqc,20)
subplot(2,3,3)
hist(psdresults{1}.smsigmasqeg,20)

subplot(2,3,4)
hist(sigmasqa_mle,20)
subplot(2,3,5)
hist(sigmasqc_mle,20)
subplot(2,3,6)
hist(sigmasqe_mle,20)


mean(sigmasqa_mle)
mean(sigmasqc_mle)
mean(sigmasqe_mle)


load('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/mwle_CortThick676_HIGH.mat')
smsigmasqa_mwle = [mwleresultsL_HIGH.sigmasqA mwleresultsR_HIGH.sigmasqA];
smsigmasqc_mwle = [mwleresultsL_HIGH.sigmasqC mwleresultsR_HIGH.sigmasqC];
smsigmasqe_mwle = [mwleresultsL_HIGH.sigmasqE mwleresultsR_HIGH.sigmasqE];
h2_mwle = [mwleresultsL_HIGH.h2' mwleresultsR_HIGH.h2'];
input = [smsigmasqa_mwle',smsigmasqc_mwle',smsigmasqe_mwle',h2_mwle'];
outfile = '~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/heritability_mwle';
myft_write_cifti(input,outfile);

mean(h2_mwle)
std(h2_mwle)


sigmag_mwle = smsigmasqa_mwle  + smsigmasqc_mwle + smsigmasqe_mwle -   sigmaemresults{1}.smsigmasqemLR;
h2_mwle_mec = smsigmasqa_mwle./sigmag_mle;
mean(h2_mwle_mec)
std(h2_mwle_mec)


%% Compare heritability estimates to sandwich truncated estimators:
% equivalent to average variance in truncated estimators:444

sum(swresults{1}.outfull_swLR.valSA(swresults{1}.outfull_swLR.valSA>0))/size(swresults{1}.outfull_swLR.vecSA,1)
sum(swresults{1}.outfull_swLR.valSC(swresults{1}.outfull_swLR.valSC>0))/size(swresults{1}.outfull_swLR.vecSC,1)
% truncation has large bias


% NOTE: Used power-iterations to estimate N eigenvalues.
sum(swresults{1}.outfull_swLR.valSA<0)
sum(swresults{1}.outfull_swLR.valSC<0)
sum(swresults{1}.outfull_swLR.valSEg<0)
sum(swresults{1}.outfull_swLR.valSEg(swresults{1}.outfull_swLR.valSEg<0)) % negligible contribution of negative eigenvalues

for v=1:nVertex
    swresults{1}.smsigmasqa(v) = swresults{1}.outfull_swLR.vecSA(v,:)*diag(swresults{1}.outfull_swLR.valSA)*swresults{1}.outfull_swLR.vecSA(v,:)';
    swresults{1}.smsigmasqc(v) = swresults{1}.outfull_swLR.vecSC(v,:)*diag(swresults{1}.outfull_swLR.valSC)*swresults{1}.outfull_swLR.vecSC(v,:)';
    swresults{1}.smsigmasqeg(v) = swresults{1}.outfull_swLR.vecSEg(v,:)*diag(swresults{1}.outfull_swLR.valSEg)*swresults{1}.outfull_swLR.vecSEg(v,:)';
end
swresults{1}.h2 = swresults{1}.smsigmasqa./(swresults{1}.smsigmasqa + swresults{1}.smsigmasqc + swresults{1}.smsigmasqeg);





mean(psdresults{1}.smsigmasqa)
mean(psdresults{1}.smsigmasqc)
mean(psdresults{1}.smsigmasqeg)

mean(psdresults{1}.h2)
std(psdresults{1}.h2)
sum(psdresults{1}.h2==0)
max(psdresults{1}.h2)
min(psdresults{1}.h2)

mean(swresults{1}.h2)
std(swresults{1}.h2)
max(swresults{1}.h2)
min(swresults{1}.h2)
 
mean(h2_mle)
std(h2_mle)

mean(h2_smmle)
std(h2_smmle)


% mean(h2_mwle)
% std(h2_mwle)

hist(h2_mwle)
median(h2_mwle)

%% Compare variance for S-FSEM:

load('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_heritabilities.mat')
h2_sfsem = [out_swL.h2_symm;out_swR.h2_symm];
mean(h2_sfsem)
min(h2_sfsem)
max(h2_sfsem)
std(h2_sfsem)

input = [[out_swL.smSA_symm_variance;out_swR.smSA_symm_variance],[out_swL.smSC_symm_variance;out_swR.smSC_symm_variance],...
    [out_swL.smSEg_symm_variance;out_swR.smSEg_symm_variance],[out_swL.h2_symm;out_swR.h2_symm]];
mean(input)




%% Compare variance from covariance function to variances from smoothed MLE
figure;
subplot(2,2,1)
scatter(smsigmasqa_smmle,psdresults{1}.smsigmasqa)
ylabel('\sigma_a PSD-ACE')
xlabel('\sigma_a MLE')
hline = refline(1);
hline.Color = 'k';
axis equal;
title('PSD-ACE versus MLE: \sigma_a')

subplot(2,2,2)
scatter(smsigmasqc_smmle,psdresults{1}.smsigmasqc)
ylabel('PSD-ACE')
hline = refline(1);
hline.Color = 'k';
axis equal;
title('PSD-ACE versus MLE: \sigma_c')

subplot(2,2,3)
scatter(smsigmasqe_smmle,psdresults{1}.smsigmasqeg)
ylabel('PSD-ACE')
hline = refline(1);
hline.Color = 'k';
axis equal;
title('PSD-ACE versus MLE: \sigma_{e,G}')


psdace_total = psdresults{1}.smsigmasqa + psdresults{1}.smsigmasqc + psdresults{1}.smsigmasqeg + sigmaemresults{1}.smsigmasqemLR';
smle_total = smsigmasqa_smmle+smsigmasqc_smmle+smsigmasqe_smmle;
subplot(2,2,4)
scatter(smle_total,psdace_total)
ylabel('PSD-ACE')
hline = refline(1);
hline.Color = 'k';
axis equal;
title('PSD-ACE versus MLE: \sigma_{T}');

saveas(gca,'~/Dropbox/SpatialACDE/Documents/Figures/ScatterPSDACEvSMLE.jpg')

mean(psdace_total)
mean(smle_total)



%% Create subsets covariance matrices

% for input to myft_write_seedcorr:
psdresults{1}.smXa = psdresults{1}.outfull_psdLR.Xa;

myvertex = 32512+1; %reminder: wb_view numbers from 0
tempindex = find(mapping_CORTEX_RIGHT==(myvertex-32492));
myindices = mapping_CORTEX_RIGHT(tempindex:tempindex+499)+32492;

% write cifti file:
myseedcorr = myft_write_seedcorr(myindices,psdresults,...
    ['~/Dropbox/SpatialACDE/Results/CorrSigmaA_conpsd_seeds',num2str(min(myindices)),'to',num2str(max(myindices))]);

% write cifti file:
myft_write_seedcorr(30402:30501,psdresults,'~/Dropbox/SpatialACDE/Results/CorrSigmaA_conpsd_seeds30402to30501',...
    '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii',...
    '~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L',...
    '~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R');



