%----------------------
% Examine point-wise MLE for cortical thickness
% NOTE: This file conducts variable selection and
% creates diagnostic figures and is intended to be run
% locally.
% Updated 10 April 2017: Re-ran with 1200-subject data release
%------------------------------------

load '~/Dropbox/SpatialACDE/Data/subjectData_cortThickForMLE.mat';
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R
load ~/Dropbox/SpatialACDE/Data/surfdataLR.mat

addpath '~/Dropbox/SpatialACDE/Programs/Functions'
nSubject = size(dataMatL,2);

%% Calculate Vertex-wise MLE:
% also conduct variable selection
%   variable selection conducted in 5_EstimateACE_cortThick_LOW.m

tic;
mleresultsR_HIGH = allvertexmle_ACE(dataMatR', xMat, familyst);
toc
% note: the warnings of matrices close to singular are from the
% inversion of the fisher information matrix. This occurs
% when the variance components are estimated to be equal to zero,
% as they are included in the FI. The approximate p-values from the FI
% for the fixed effects appear unaffected. 

covnames = {'intercept','gender','age','handedness','height','weight','BMI','ICV'};

figure;
for i=1:8
    subplot(3,3,i)
    hist(mleresultsR_HIGH.pvalueBeta(i,:))
    title(covnames(i));
end
saveas(gca,'~/Dropbox/SpatialACDE/Figures/Examine_pvalues_allcovariates_cortThick.jpg')
save('~/Dropbox/SpatialACDE/Results/mle_cortThick676_SelectCovariates.mat','mleresultsR_HIGH')



tabulate(mleresultsR_HIGH.eFlag_a)
% values should be greater than 0; 0 means maxit reached before convergence
tabulate(mleresultsR_HIGH.eFlag_n)

figure;
for i=1:4
    subplot(2,2,i)
    hist(mleresultsR_HIGH.pvalueBeta(i,:))
    title(selcovnames(i));
end



%% Not run -------->
if 0
%% Create diagnostic plots:
rng(321)
subIndex = sort(datasample(1:nVertexR,25));
a = figure;
for k=1:25
    subplot(5,5,k)
    scatter(mleresultsR_HIGH.yfixed(:,subIndex(k)),mleresultsR_HIGH.fixefresid(:,subIndex(k)),'.')
    hline = refline(0,0);
    hline.Color = 'r';
    title(['Vertex: ',num2str(subIndex(k))],'Fontsize',9)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 6)
end
suptitle('Residuals vs Fitted, right hemisphere')
saveas(a,'~/Dropbox/SpatialACDE/Figures/ResidsVFitted_CortThickR.png');

% TO DO: Sometimes the fitted values have a bimodal distribution due
% to a large gender effect, but this is not apparent in the raw data.
% why is the effect size sometimes overestimated.

for k=1:25
    subplot(5,5,k)
    hist(mleresultsR_HIGH.yfixed(:,subIndex(k)))
    title(['Vertex: ',num2str(k)],'Fontsize',9)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 6)
end
suptitle('Histograms, right hemisphere')
saveas(a,'~/Dropbox/SpatialACDE/Figures/Histograms_CortThickR.png');


hist(mleresultsR_HIGH.pvalueBeta(2,:))
hist(mleresultsR_HIGH.pvalueBeta(3,:))

% write results to cifti:
% resultsThickL= [mleresultsL_HIGH.heritability', mleresultsL_HIGH.sigmasqA', mleresultsL_HIGH.sigmasqC', mleresultsL_HIGH.sigmasqE'];
% createFuncCIItoGIIgeneric(resultsThickL,'CORTEX_LEFT','~/Dropbox/SpatialACDE/Results/cortThickLVarianceEstimates556_mle')
% 
% resultsThickR= [mleresultsR_HIGH.heritability', mleresultsR_HIGH.sigmasqA', mleresultsR_HIGH.sigmasqC', mleresultsR_HIGH.sigmasqE'];
% createFuncCIItoGIIgeneric(resultsThickR,'CORTEX_RIGHT','~/Dropbox/SpatialACDE/Results/cortThickRVarianceEstimates556_mle')

%----------------------------






%%% plot surfaces in matlab

gii_smmleR_HIGH = zeros(32492,1);
gii_smmleR_HIGH(mapping_CORTEX_RIGHT) = smmleR_HIGH.smsigmasqa;


gii_varAfullR = zeros(32492,1);
gii_varAfullR(mapping_CORTEX_RIGHT) = mleresultsR_HIGH.sigmasqA;


sphere_R = gifti('~/Dropbox/MyHCP/Data/100307/MNINonLinear/fsaverage_LR32k/100307.R.sphere.32k_fs_LR.surf.gii');


h=subplot(2,2,1);
trisurf(surfdata_R.faces, surfdata_R.vertices(:,1), surfdata_R.vertices(:,2), surfdata_R.vertices(:,3),...
        gii_varAfullR,'EdgeColor','None');
axis equal
axis off
colormap('jet')
view(-90,0);
colorrange = caxis();
caxis([0 0.1])
    
h=subplot(2,2,2);
trisurf(surfdata_R.faces, surfdata_R.vertices(:,1), surfdata_R.vertices(:,2), surfdata_R.vertices(:,3),...
        gii_varAfullR,'EdgeColor','None');
axis equal
axis off
colormap('jet')
view(90,0);
colorrange = caxis();
caxis([0 0.1])
    
h=subplot(2,2,3);
trisurf(surfdata_R.faces, surfdata_R.vertices(:,1), surfdata_R.vertices(:,2), surfdata_R.vertices(:,3),...
        gii_smmleload '~/Dropbox/SpatialACDE/Results/smmle_cortThick556_HIGH.mat'
R_HIGH,'EdgeColor','None');
axis equal
axis off
colormap('jet')
view(-90,0);
colorrange = caxis();
caxis([0 0.1])
    
h=subplot(2,2,4);
trisurf(surfdata_R.faces, surfdata_R.vertices(:,1), surfdata_R.vertices(:,2), surfdata_R.vertices(:,3),...
        gii_smmleR_HIGH,'EdgeColor','None');
axis equal
axis off
colormap('jet')
view(90,0);
colorrange = caxis();
caxis([0 0.1])
end   
    



