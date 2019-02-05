%----------------
% Benjamin Risk
% Create surface plots of covariance features
%
 
addpath './Functions'
%addpath '~/Dropbox/SpatioTemporalFMRI/Programs/mfunctions'
% To create tables, download this matlab package:
addpath ~/Dropbox/mfunctions/eliduenisch-latexTable-5212622

% covariance matrices used in simulations:
load './Simulations/SupportingDataFiles/CovarianceMatricesForSimulationsSharmACEM_compress.mat'
% expand covariance matrices
sigmaa = sigmaa_v*sigmaa_l*sigmaa_v';
sigmac = sigmac_v*sigmac_l*sigmac_v';
sigmaeg = sigmaeg_v*sigmaeg_l*sigmaeg_v';


% simulation results:
fprintf('EDIT THIS------------->')
fprintf('load ~/Dropbox/SpatialACDE/Results/aSimulations_100_100_200_08-<date>.mat')
fprintf('<-----------')

load '~/Dropbox/SpatialACDE/Results/aSimulations_100_100_200_08-<date>.mat'
load ~/Dropbox/SpatialACDE/Results/aSimulations_100_100_200_08-Feb-2017_em0pt2.mat


%depthLabels = {'Sandwich Symm','Sandwich Trunc','SL Symm','SL Trunc','PSD Est Rank','PSD True Rank','MLE','SMLE','MWLE'};
depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};

tabulate(allconvergence_estrank)


tabulate(allconvergence_truerank)

% true rank tends to have more flag==-2
gradnorm_truerank(allconvergence_truerank==-2)
norm(sigmaa,'fro')
% these are acceptable


% summary measures that appear in the manuscript:
% PSD-ACE estimates are less variable:
var(diag(sigmaa))
var(diag(squeeze(simresultsCovfunAVE(1,5,:,:))))
var(diag(squeeze(simresultsCovfunAVE(1,6,:,:))))

mean(diag(sigmaa))

mean(diag(squeeze(simresultsCovfunAVE(1,5,:,:)))) % PSD-ACE
mean(diag(squeeze(simresultsCovfunAVE(1,6,:,:)))) % PSD-ACE-O
% biased upwards 

mean(diag(squeeze(simresultsCovfunAVE(1,2,:,:)))) % PSD-SW
mean(diag(squeeze(simresultsCovfunAVE(1,4,:,:)))) %PSD-FSEM
% biased upwards 


mean(diag(squeeze(simresultsCovfunAVE(1,1,:,:)))) % S-SW
mean(diag(squeeze(simresultsCovfunAVE(1,3,:,:)))) %S-FSEM
% little bias

mean(diag(squeeze(simresultsCovfunAVE(1,7,:,:)))) %MLE
mean(diag(squeeze(simresultsCovfunAVE(1,9,:,:)))) %MWLE
mean(diag(squeeze(simresultsCovfunAVE(1,8,:,:)))) %SMLE
% little bias

%--------sigmasqC
mean(diag(squeeze(simresultsCovfunAVE(2,5,:,:)))) % PSD-ACE
mean(diag(squeeze(simresultsCovfunAVE(2,6,:,:)))) % PSD-ACE-O
% biased upwards 

mean(diag(squeeze(simresultsCovfunAVE(2,2,:,:)))) % PSD-SW
mean(diag(squeeze(simresultsCovfunAVE(2,4,:,:)))) %PSD-FSEM
% biased upwards 


mean(diag(squeeze(simresultsCovfunAVE(2,1,:,:)))) % S-SW
mean(diag(squeeze(simresultsCovfunAVE(2,3,:,:)))) %S-FSEM
% unbiased

mean(diag(squeeze(simresultsCovfunAVE(2,7,:,:)))) %MLE
mean(diag(squeeze(simresultsCovfunAVE(2,9,:,:)))) %MWLE
mean(diag(squeeze(simresultsCovfunAVE(2,8,:,:)))) %SMLE
% biased upwards 




%%-------------------------------------------------
% summarize Bias^2, Variance, and MSE of covariance functions:
%
% NOTE: This results appear in 
input.dataFormat = {'%.2f'}; 
reSortIndexSix = [3,4,1,2,6,5];
input.data = [simresultsCovMatBiassq_sumV(1,reSortIndexSix)',simresultsCovMatVariance_sumV(1,reSortIndexSix)',...
    simresultsCovMatMSE_sumV(1,reSortIndexSix)',simresultsCovMatBiassq_sumV(2,reSortIndexSix)',...
    simresultsCovMatVariance_sumV(2,reSortIndexSix)',simresultsCovMatMSE_sumV(2,reSortIndexSix)',...
    simresultsCovMatBiassq_sumV(3,reSortIndexSix)',simresultsCovMatVariance_sumV(3,reSortIndexSix)',...
    simresultsCovMatMSE_sumV(3,reSortIndexSix)'];

input.tableRowLabels = depthLabels(reSortIndexSix);
latexTable(input);

%%-----------------------+------------------------
% summarize Bias^2, Variance, and MSE for variance parameters:
reSortIndexNine = [7,9,8,3,4,1,2,6,5];
input.dataFormat = {'%.2f'}; 
input.data = [simresultsVarVecBiassq_sumV(1,reSortIndexNine)',simresultsVarVecVariance_sumV(1,reSortIndexNine)',...
    simresultsVarVecMSE_sumV(1,reSortIndexNine)',simresultsVarVecBiassq_sumV(2,reSortIndexNine)',...
    simresultsVarVecVariance_sumV(2,reSortIndexNine)',simresultsVarVecMSE_sumV(2,reSortIndexNine)',...
    simresultsVarVecBiassq_sumV(3,reSortIndexNine)',simresultsVarVecVariance_sumV(3,reSortIndexNine)',...
    simresultsVarVecMSE_sumV(3,reSortIndexNine)',...
    simresultsh2Biassq_sumV(reSortIndexNine),simresultsh2Variance_sumV(reSortIndexNine),simresultsh2MSE_sumV(reSortIndexNine),];

%input.transposeTable = 1;
input.tableRowLabels = depthLabels(reSortIndexNine);
latexTable(input);





%-----------------------------
% use the same color scale for all images:

colorrange = [0.01,0.3];

figure(1);
c=subplot(4,3,1);
trisurf(face,x,y,z,h2,'EdgeColor','none')
%colormap('jet')
%colorbar()
title('Truth')
%colorrange=caxis();
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,2)
h=trisurf(face,x,y,z,h2,'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;

c=subplot(4,3,4);
trisurf(face,x,y,z,simresultsh2AVE(7,:),'EdgeColor','none')
%colorbar()
title(depthLabels{7})
caxis(colorrange);
axis equal;
axis off;

c=subplot(4,3,5);
trisurf(face,x,y,z,simresultsh2AVE(9,:),'EdgeColor','none')
%colorbar()
title(depthLabels{9})
caxis(colorrange)
axis equal;
axis off;

c=subplot(4,3,6);
trisurf(face,x,y,z,simresultsh2AVE(8,:),'EdgeColor','none')
%colorbar()
title(depthLabels{8})
caxis(colorrange);
axis equal;
axis off;

subplot(4,3,7)
trisurf(face,x,y,z,simresultsh2AVE(3,:),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{3})
axis equal;
axis off;

subplot(4,3,8)
trisurf(face,x,y,z,simresultsh2AVE(4,:),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title(depthLabels{4})
axis equal;
axis off;

subplot(4,3,9)
trisurf(face,x,y,z,simresultsh2AVE(1,:),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title(depthLabels{1})
axis equal;
axis off;

subplot(4,3,10)
trisurf(face,x,y,z,simresultsh2AVE(2,:),'EdgeColor','none')
caxis(colorrange)
title(depthLabels{2})
axis equal;
axis off;

subplot(4,3,11)
trisurf(face,x,y,z,simresultsh2AVE(6,:),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{6})
axis equal;
axis off;

subplot(4,3,12)
trisurf(face,x,y,z,simresultsh2AVE(5,:),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{5})
axis equal;
axis off;
saveas(c,'~/Dropbox/SpatialACDE/Documents/Figures/Simulations_h2_100_100_200_fixedcaxis.png')

%hp4 = get(subplot(3,3,3),'Position')
%colorbar('Position', [0.7,7.01,0.2,0.201])

%%----------------------------------------------------
%-----------------------------
% use matlab-chosen color scale for each figure:
figure(2);
c=subplot(4,3,1);
trisurf(face,x,y,z,h2,'EdgeColor','none')
colorbar()
title('Truth')
axis equal;
axis off;

subplot(4,3,4);
trisurf(face,x,y,z,simresultsh2AVE(7,:),'EdgeColor','none')
colorbar()
title(depthLabels{7})
axis equal;
axis off;

subplot(4,3,5);
trisurf(face,x,y,z,simresultsh2AVE(9,:),'EdgeColor','none')
colorbar()
title(depthLabels{9})
axis equal;
axis off;

subplot(4,3,6);
trisurf(face,x,y,z,simresultsh2AVE(8,:),'EdgeColor','none')
colorbar()
title(depthLabels{8})
axis equal;
axis off;

subplot(4,3,7)
trisurf(face,x,y,z,simresultsh2AVE(3,:),'EdgeColor','none')
colorbar()
title(depthLabels{3})
axis equal;
axis off;

subplot(4,3,8)
trisurf(face,x,y,z,simresultsh2AVE(4,:),'EdgeColor','none')
colorbar()
title(depthLabels{4})
axis equal;
axis off;

subplot(4,3,9)
trisurf(face,x,y,z,simresultsh2AVE(1,:),'EdgeColor','none')
colorbar()
title(depthLabels{1})
axis equal;
axis off;

subplot(4,3,10)
trisurf(face,x,y,z,simresultsh2AVE(2,:),'EdgeColor','none')
colorbar()
title(depthLabels{2})
axis equal;
axis off;

subplot(4,3,11)
trisurf(face,x,y,z,simresultsh2AVE(6,:),'EdgeColor','none')
 
colorbar()
title(depthLabels{6})
axis equal;
axis off;

subplot(4,3,12)
trisurf(face,x,y,z,simresultsh2AVE(5,:),'EdgeColor','none')
 
colorbar()
title(depthLabels{5})
axis equal;
axis off;
saveas(c,'~/Dropbox/SpatialACDE/Documents/Figures/Simulations_h2_100_100_200_varyingcaxis.png')


%--------------------

% test: look at images for average of sigma a / (average of sigma a +
% average of sigma c + average of sigma e)

a = diag(squeeze(simresultsCovfunAVE(1,5,:,:)));
c = diag(squeeze(simresultsCovfunAVE(2,5,:,:)));
e = simresultsSigmasqeAVE(2,:)';
alt_h2ave = a./(a+c+e);


subplot(2,1,1)
trisurf(face,x,y,z,simresultsh2AVE(5,:),'EdgeColor','none')
 
colorbar()
title(depthLabels{5})
axis equal;
axis off;

subplot(2,1,2)
trisurf(face,x,y,z,alt_h2ave,'EdgeColor','none')
 
colorbar()
title(depthLabels{5})
axis equal;
axis off;


%-----------------------------------------------
%----------------------------------
%% Sigma A: ---------------------------

c=figure(3);
subplot(4,3,1);
trisurf(face,x,y,z,diag(sigmaa),'EdgeColor','none')
 
%colorbar()
title('Truth')
%colorrange=caxis();
colorrange = [0,0.0375];
caxis(colorrange);
axis equal;
axis off;

subplot(4,3,2)
h=trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,6,:,:))),'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;

subplot(4,3,4);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,7,:,:))),'EdgeColor','none')
 
%colorbar()
title('MLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,5);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,9,:,:))),'EdgeColor','none')
 
%colorbar()
title('MWLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,6);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,8,:,:))),'EdgeColor','none')
 
%colorbar()
title('SMMLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,7)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,3,:,:))),'EdgeColor','none')
 
caxis(colorrange);
title('S-FSEM')
axis equal;
axis off;

subplot(4,3,8)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,4,:,:))),'EdgeColor','none')
 
%colorbar()
caxis(colorrange);
title('PSD-FSEM')
axis equal;
axis off;

subplot(4,3,9)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,1,:,:))),'EdgeColor','none')
 
%colorbar()
caxis(colorrange);
title('S-SW')
axis equal;
axis off;

subplot(4,3,10)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,2,:,:))),'EdgeColor','none')
 
caxis(colorrange)
title('PSD-SW')
axis equal;
axis off;

subplot(4,3,11)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,5,:,:))),'EdgeColor','none')
 
caxis(colorrange);
title('PSD-ACE-O')
axis equal;
axis off;

subplot(4,3,12)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(1,6,:,:))),'EdgeColor','none')
 
caxis(colorrange);
title('PSD-ACE')
axis equal;
axis off;

%hp4 = get(subplot(3,3,3),'Position')
%colorbar('Position', [0.7,7.01,0.2,0.201])

saveas(c,'~/Dropbox/SpatialACDE/Documents/Figures/Simulations_SigmasqA_100_100_200.png')

%% -------------------------------------------
%% Sigma C
figure(4);
d=subplot(4,3,1);
trisurf(face,x,y,z,diag(sigmac),'EdgeColor','none')
%colormap('jet')
%colorbar()
title('Truth')
colorrange = [0,0.03];
%colorrange=caxis();
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,2)
h=trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,6,:,:))),'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;

subplot(4,3,4);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,7,:,:))),'EdgeColor','none')
%colorbar()
title('MLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,5);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,9,:,:))),'EdgeColor','none')
%colorbar()
title('MWLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,6);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,8,:,:))),'EdgeColor','none')
%colorbar()
title('SMMLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,7)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,3,:,:))),'EdgeColor','none')
caxis(colorrange);
title('S-FSEM')
axis equal;
axis off;

subplot(4,3,8)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,4,:,:))),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title('PSD-FSEM')
axis equal;
axis off;

subplot(4,3,9)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,1,:,:))),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title('S-SW')
axis equal;
axis off;

subplot(4,3,10)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,2,:,:))),'EdgeColor','none')
caxis(colorrange)
title('PSD-SW')
axis equal;
axis off;

subplot(4,3,11)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,5,:,:))),'EdgeColor','none')
caxis(colorrange);
title('PSD-ACE-O')
axis equal;
axis off;

subplot(4,3,12)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(2,6,:,:))),'EdgeColor','none')
caxis(colorrange);
title('PSD-ACE')
axis equal;
axis off;

saveas(d,'~/Dropbox/SpatialACDE/Documents/Figures/Simulations_SigmasqC_100_100_200.png')

%% -------------------------------------------
%% Sigma Eg
figure(5);
d=subplot(4,3,1);
trisurf(face,x,y,z,diag(sigmaeg),'EdgeColor','none')
%colormap('jet')
%colorbar()
title('Truth')
%colorrange = [0,0.03];
colorrange=caxis();
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,2)
h=trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,6,:,:))),'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;

subplot(4,3,4);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,7,:,:))),'EdgeColor','none')
%colorbar()
title('MLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,5);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,9,:,:))),'EdgeColor','none')
%colorbar()
title('MWLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,6);
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,8,:,:))),'EdgeColor','none')
%colorbar()
title('SMMLE')
caxis(colorrange)
axis equal;
axis off;

subplot(4,3,7)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,3,:,:))),'EdgeColor','none')
caxis(colorrange);
title('S-FSEM')
axis equal;
axis off;

subplot(4,3,8)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,4,:,:))),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title('PSD-FSEM')
axis equal;
axis off;

subplot(4,3,9)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,1,:,:))),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title('S-SW')
axis equal;
axis off;

subplot(4,3,10)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,2,:,:))),'EdgeColor','none')
caxis(colorrange)
title('PSD-SW')
axis equal;
axis off;

subplot(4,3,11)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,5,:,:))),'EdgeColor','none')
caxis(colorrange);
title('PSD-ACE-O')
axis equal;
axis off;

subplot(4,3,12)
trisurf(face,x,y,z,diag(squeeze(simresultsCovfunAVE(3,6,:,:))),'EdgeColor','none')
caxis(colorrange);
title('PSD-ACE')
axis equal;
axis off;

saveas(d,'~/Dropbox/SpatialACDE/Documents/Figures/Simulations_SigmasqEg_100_100_200.png')





%%----------------------------------------------------
%% take a row of the covariance matrix of Sigma A:
rng(321)

iIndex = randsample(1002,1);


figure(3);
h=subplot(3,3,1);
trisurf(face,x,y,z,sigmaa(:,iIndex),'EdgeColor','none')
%colorbar()
title('Truth')
colorrange=caxis();
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
mypos =  [-22.9690   -0.2004   19.2740];
h.CameraPosition = mypos;

subplot(3,3,2)
h=trisurf(face,x,y,z,sigmaa(:,iIndex),'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;
axis equal;


h=subplot(3,3,4);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(1,3,:,iIndex)),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title(depthLabels{3})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;



h=subplot(3,3,5);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(1,4,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{4})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;



h=subplot(3,3,6);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(1,1,:,iIndex)),'EdgeColor','none')
caxis(colorrange)
title(depthLabels{1})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,7);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(1,2,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{2})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;



h=subplot(3,3,8);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(1,6,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{6})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,9);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(1,5,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{5})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;

saveas(h,['~/Dropbox/SpatialACDE/Documents/Figures/Simulations_SigmaA_v=' num2str(iIndex) '.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% take a row of the covariance matrix SigmaC:

%iIndex = randsample(1002,1);

figure(4);
h=subplot(3,3,1);
trisurf(face,x,y,z,sigmac(:,iIndex),'EdgeColor','none')
%colorbar()
title('Truth')
colorrange=caxis();
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
%mypos =  [-22.9240  -16.9629    9.2659];
mypos =  [-22.9690   -0.2004   19.2740];
h.CameraPosition = mypos;

subplot(3,3,2)
h=trisurf(face,x,y,z,sigmac(:,iIndex),'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;
axis equal;


h=subplot(3,3,4);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(2,3,:,iIndex)),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title(depthLabels{3})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;



h=subplot(3,3,5);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(2,4,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{4})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;



h=subplot(3,3,6);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(2,1,:,iIndex)),'EdgeColor','none')
caxis(colorrange)
title(depthLabels{1})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,7);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(2,2,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{2})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;



h=subplot(3,3,8);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(2,6,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{6})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;




h=subplot(3,3,9);
trisurf(face,x,y,z,squeeze(simresultsCovfunAVE(2,5,:,iIndex)),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{5})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;

saveas(h,['~/Dropbox/SpatialACDE/Documents/Figures/Simulations_SigmaC_v=' num2str(iIndex) '.png'])
