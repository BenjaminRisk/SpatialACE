% Estimate covariance functions using Luo Symm

addpath '~/Dropbox/SpatialACDE/Programs/Functions'
addpath ~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8

maxNumCompThreads(24);

load ~/Dropbox/SpatialACDE/Data/subjectData_cortThickForCovfun.mat;
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat


xMatSub = xMat(:,[1 2 3 8]);
selcovnames = {'intercept','gender','age','ICV'};
 
residL = dataMatL' - xMatSub*smmleL_HIGH.smbeta;
residR = dataMatR' - xMatSub*smmleR_HIGH.smbeta;
nVertexL = length(latL);
nVertexR = length(latR);



out_swL = fullcovacem_sl_symm_bigdata(residL,familyst,latL,longL,1.3);
out_swR = fullcovacem_sl_symm_bigdata(residR,familyst,latR,longR,1.3);

%out_swL = fullcovacem_sl_symm_bigdata(residL(:,1:1000),familyst,latL(1:1000),longL(1:1000),10);
%out_swR = fullcovacem_sl_symm_bigdata(residR(:,1:1000),familyst,latR(1:1000),longR(1:1000),10);
% sum([out_swL.h2_symm;out_swR.h2_symm]<0)

save('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_heritabilities.mat','out_swL','out_swR');

input = [[out_swL.smSA_symm_variance;out_swR.smSA_symm_variance],[out_swL.smSC_symm_variance;out_swR.smSC_symm_variance],...
    [out_swL.smSEg_symm_variance;out_swR.smSEg_symm_variance],[out_swL.h2_symm;out_swR.h2_symm]];

outfile = '~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/heritability_SFSEM'; 
myft_write_cifti(input,outfile);

h2 = [out_swL.h2_symm;out_swR.h2_symm];
mean(h2)
min(h2)
max(h2)
std(h2)

