%-------------------
% Bind subject cortex data and create
% vectors from the .csv output from a2_CreateSubjectListFamilyID.R
%
% NOTE: The subjects in subjList.csv have been ordered by 
% numeric zygosity (=1 if MZ, 2 if DZ, and 3 if singleton),
% then by familyID. This script organizes the subjects in the
% main data matrix in this order. This ordering is essential
% for the correct calculation of the likelihood that occurs 
% when using subsequent programs.
%
% Updates
%   - create two sets of data:
%       1. dataset for maximum likelihood excluding siblings
%       2. dataset for covfun including all subjects, 
%           siblings labeled as MDti
%---------------------------------------------
addpath ~/Dropbox/SpatialACDE/Programs/Functions

subjCov = readtable('~/Dropbox/SpatialACDE/Data/subjCovariatesForMLE.csv');
subjCovNumeric = table2array(subjCov);
nSubject = size(subjCovNumeric,1);
covnames = {'intercept','gender','age','handedness','height','weight','BMI','ICV'};

% design matrix with intracranial volume:
xMat = [ones(nSubject,1),subjCovNumeric(:,5:11)];

subjList = readtable('~/Dropbox/SpatialACDE/Data/subjRelatednessForMLE.csv');
subjListNumeric = table2array(subjList);

%create vectors for all relatedness columns:
subjID = subjListNumeric(:,1);
familyID = subjListNumeric(:,2);
familyst.MZtp1 = logical(subjListNumeric(:,3));
familyst.MZtp2 = logical(subjListNumeric(:,4));
familyst.DZtp1 = logical(subjListNumeric(:,5));
familyst.DZtp2 = logical(subjListNumeric(:,6));
familyst.MDti = logical(subjListNumeric(:,7));
familyst.familyID = familyID;

load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L.mat
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R.mat

nVertexL = length(mapping_CORTEX_LEFT);
nVertexR = length(mapping_CORTEX_RIGHT);

%addpath ~/Applications2/robertoostenveld-cifti-matlab-27383b8
% load data for subjects in subjID:
dataMatL = zeros(nVertexL,nSubject);
dataMatR = zeros(nVertexR,nSubject);

tic;
for n=1:nSubject
    dataMatL(:,n) = readincortexst(['~/Dropbox/MyHCP/Data/cortThick32kLR_thickness_MSMAll/' num2str(subjID(n)) '.thickness_MSMAll.32k_fs_LR.dscalar.nii'],'CORTEX_LEFT');
    dataMatR(:,n) = readincortexst(['~/Dropbox/MyHCP/Data/cortThick32kLR_thickness_MSMAll/' num2str(subjID(n)) '.thickness_MSMAll.32k_fs_LR.dscalar.nii'],'CORTEX_RIGHT');
end
toc

%clear subjList subjListNumeric subjCovNumeric nSubject nVertexL nVertexR mapping_CORTEX_LEFT mapping_CORTEX_RIGHT
save('~/Dropbox/SpatialACDE/Data/subjectData_cortThickForMLE','familyst','xMat','dataMatL','dataMatR','covnames')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same code for creating second dataset:
subjCov = readtable('~/Dropbox/SpatialACDE/Data/subjCovariatesForCovfun.csv');
subjCovNumeric = table2array(subjCov);
nSubject = size(subjCovNumeric,1);

% design matrix with intracranial volume:
xMat = [ones(nSubject,1),subjCovNumeric(:,5:11)];

subjList = readtable('~/Dropbox/SpatialACDE/Data/subjRelatednessForCovfun.csv');
subjListNumeric = table2array(subjList);

%create vectors for all relatedness columns:
subjID = subjListNumeric(:,1);
familyID = subjListNumeric(:,2);
familyst.MZtp1 = logical(subjListNumeric(:,3));
familyst.MZtp2 = logical(subjListNumeric(:,4));
familyst.DZtp1 = logical(subjListNumeric(:,5));
familyst.DZtp2 = logical(subjListNumeric(:,6));
familyst.MDti = logical(subjListNumeric(:,7));
familyst.familyID = familyID;

dataMatL = zeros(nVertexL,nSubject);
dataMatR = zeros(nVertexR,nSubject);

tic;
for n=1:nSubject
    dataMatL(:,n) = readincortexst(['~/Dropbox/MyHCP/Data/cortThick32kLR_thickness_MSMAll/' num2str(subjID(n)) '.thickness_MSMAll.32k_fs_LR.dscalar.nii'],'CORTEX_LEFT');
    dataMatR(:,n) = readincortexst(['~/Dropbox/MyHCP/Data/cortThick32kLR_thickness_MSMAll/' num2str(subjID(n)) '.thickness_MSMAll.32k_fs_LR.dscalar.nii'],'CORTEX_RIGHT');
end
toc

save('~/Dropbox/SpatialACDE/Data/subjectData_cortThickForCovfun','familyst','xMat','dataMatL','dataMatR','covnames')

