%----------------------
% Performs steps 1 to 4 in spatial ACE estimation
%
%------------------------------------
addpath '~/Dropbox/SpatialACDE/Programs/Functions'


estStep1and2 = 1;
estStep3 = 1;


%% Step 1:
if estStep1and2
    
fprintf('\nBeginning Step 1\n')

tic;

load '~/Dropbox/SpatialACDE/Data/subjectData_cortThickForMLE.mat';
% dataset for MLE: contains MZs, DZs, and randomly selected singletons
% dataset for cov: all subjects.
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat


xMatSub = xMat(:,[1 2 3 8]);
selcovnames = {'intercept','gender','age','ICV'};

parallel.importProfile('~/Dropbox/SpatialACDE/Programs/local_Copy.settings')
poolobj = parpool('local_Copy',48);
poolobj.NumWorkers

%matlabpool open 64

mleresultsR_HIGH = allvertexmle_ACE(dataMatR', xMatSub, familyst);
mleresultsL_HIGH = allvertexmle_ACE(dataMatL', xMatSub, familyst);

mytime = toc

mleresultsADE_R_HIGH = allvertexmle_ADE(dataMatR',xMatSub,familyst);
mleresultsADE_L_HIGH = allvertexmle_ADE(dataMatL',xMatSub,familyst);

%matlabpool close


save('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/mle_cortThick676_HIGH.mat','mleresultsR_HIGH',...
    'mleresultsL_HIGH','mleresultsADE_R_HIGH','mleresultsADE_L_HIGH','mytime');

mean(mleresultsADE_L_HIGH.n2loglik<mleresultsL_HIGH.n2loglik)
mean(mleresultsADE_R_HIGH.n2loglik<mleresultsR_HIGH.n2loglik)
    
    
if 0

    diff = mleresultsADE_L_HIGH.n2loglik-mleresultsL_HIGH.n2loglik;
    hist(diff(randsample(29000,1000)))

end
delete(poolobj)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Step 2: smmle

fprintf('\nBeginning Step 2\n')

%% Use GCV to smooth estimates:
maxNumCompThreads(24);

myhvec = [1.18,1.2:0.1:2];
tic;
smmleR_HIGH = gcv_smoothmle(mleresultsR_HIGH,familyst,latR,longR,myhvec,false,false);
smmleL_HIGH = gcv_smoothmle(mleresultsL_HIGH,familyst,latL,longL,myhvec,false,false);
toc


save('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat','smmleR_HIGH','smmleL_HIGH');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:
if estStep3
fprintf('\nBeginning Step 3\n')

% clear old datasets with 676 subjects:
clear

% Load all subjects:
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load ~/Dropbox/SpatialACDE/Data/subjectData_cortThickForCovfun.mat;

K=1;
myhvec = [1.18,1.2:0.1:2];

tic;

% calculate residuals for all 1094 subjects:
xMatSub = xMat(:,[1 2 3 8]);

residL = dataMatL' - xMatSub*smmleL_HIGH.smbeta;
residR = dataMatR' - xMatSub*smmleR_HIGH.smbeta;

smsigmasqaL = smmleL_HIGH.smsigmasqa;
smsigmasqcL = smmleL_HIGH.smsigmasqc;
smsigmasqeL = smmleL_HIGH.smsigmasqe;
smsigmasqaR = smmleR_HIGH.smsigmasqa;
smsigmasqcR = smmleR_HIGH.smsigmasqc;
smsigmasqeR = smmleR_HIGH.smsigmasqe;

nVertexL = length(latL);
nVertexR = length(latR);
rng(111)

sigmaemresults = cell(K,1);

if K>1
    cvpartL = cvpartition(nVertexL,'Kfold',K);
    cvpartR = cvpartition(nVertexR, 'Kfold',K);
else
    cvpartL = 'all';
    cvpartR = 'all';
end
    
    for k=1:K
        if K>1
            randLatR = latR(cvpartR.test(k));
            randLongR = longR(cvpartR.test(k));
            vertexSubsetR = find(cvpartR.test(k)); % get indices instead of logical vector

            randLatL = latL(cvpartL.test(k));
            randLongL = longL(cvpartL.test(k));
            vertexSubsetL = find(cvpartL.test(k));
        else
            randLatR = latR;
            randLongR = longR;
            vertexSubsetR = [1:nVertexR]';

            randLatL = latL;
            randLongL = longL;
            vertexSubsetL = [1:nVertexL]';
        end
        nR = length(vertexSubsetR);
        nL = length(vertexSubsetL);
    
        subResidL = residL(:,vertexSubsetL);
        subResidR = residR(:,vertexSubsetR);
    
        sigmaemresults{k}.estsmsigmasqemL = estsigmasqem_gcv(subResidL,familyst,smsigmasqaL(vertexSubsetL),...
            smsigmasqcL(vertexSubsetL),smsigmasqeL(vertexSubsetL),randLatL,randLongL,myhvec);
        % note: if mse = [NaN,NaN], chooses first bw; if mse = [NaN,any
        % number], chooses second bw.


        % compare with cv:
        % this will take a massive amount of time... 
        %tic;
        %sigmaemresults.estsmsigmasqemL_CV = estsigmasqem_cv(residL,smsigmasqaL,smsigmasqcL,smsigmasqeL,latL,longL,myhvec,5);
        %toc;

        fprintf('\nCompleted left hemisphere\n')

        %NOTE: The results estimating each hemisphere separately are equivalent
        %to estimating them simultaneously; this is easily seen by using block
        %matrix notation, and seeing block diagonals are K_L R_L R_L^T K_L  and
        %K_R R_R R_R^T K_R, and note off block is K_R R_R R_L^T K_L
        sigmaemresults{k}.estsmsigmasqemR = estsigmasqem_gcv(subResidR,familyst,smsigmasqaR(vertexSubsetR),...
            smsigmasqcR(vertexSubsetR),smsigmasqeR(vertexSubsetR),randLatR,randLongR,sigmaemresults{k}.estsmsigmasqemL.hvecmin);

        sigmaemresults{k}.smsigmasqemLR = [sigmaemresults{k}.estsmsigmasqemL.sigmasqem,sigmaemresults{k}.estsmsigmasqemR.sigmasqem];

        nmissing = sum(isnan(sigmaemresults{k}.smsigmasqemLR));
        if nmissing
            warning([num2str(nmissing) ' values imputed due to missingness'])
            sigmaemresults{k}.smsigmasqemLR = impute_missing_sigmasqem(sigmaemresults{k}.smsigmasqemLR,[randLatL;randLatR],...
                [randLongL;randLongR]);
        end
    end   

mytime = toc;
save('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/SigmaemCortThickACEM.mat','sigmaemresults','mytime','cvpartR','cvpartL');
%delete(poolobj);
end


%%%%%%%%%%%%%%%%%%%%%
%% Step 4: Generate preliminary estimators using SW

fprintf('\nBeginning Step 4\n')

K=1; % NOTE: K=1 uses the full dataset. K>1 partitions the data to 
	% decrease memory requirements, but is less accurate
swresults = cell(K,1);

%NOTE: Set of bandwidths need to change depending on the resolution
myhvec = [1.18,1.25,1.3,1.35,1.4,1.5,1.75];

% comment out if running full script:
if ~estStep3
    load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat
    load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
    load ~/Dropbox/SpatialACDE/Data/subjectData_cortThickForCovfun.mat;
    load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/SigmaemCortThickACEM.mat;

    % calculate residuals for all 1094 subjects:
    xMatSub = xMat(:,[1 2 3 8]);
    residL = dataMatL' - xMatSub*smmleL_HIGH.smbeta;
    residR = dataMatR' - xMatSub*smmleR_HIGH.smbeta;
    nVertexL = length(latL);
    nVertexR = length(latR);
    rng(111)
    if K>1
        cvpartL = cvpartition(nVertexL,'Kfold',K);
        cvpartR = cvpartition(nVertexR, 'Kfold',K);
    else
        cvpartL = 'all';
        cvpartR = 'all';
    end


    % comment out if performing gcv on sw:
    %myhvec = sigmaemresults{1}.estsmsigmasqemL.hvecmin;
end


nSubject = size(dataMatL,2);
neigsSA = sum(familyst.MZtp1)+sum(familyst.DZtp1)+1;
neigsSC = neigsSA;
neigsSEg = nSubject -sum(familyst.MZtp1)+10;


tic;
%% Loop through data partitions
for k=1:K
    
    if K>1
        randLatR = latR(cvpartR.test(k));
        randLongR = longR(cvpartR.test(k));
        vertexSubsetR = find(cvpartR.test(k)); % get indices instead of logical vector
        
        randLatL = latL(cvpartL.test(k));
        randLongL = longL(cvpartL.test(k));
        vertexSubsetL = find(cvpartL.test(k));
    else
        randLatR = latR;
        randLongR = longR;
        vertexSubsetR = [1:nVertexR]';

        randLatL = latL;
        randLongL = longL;
        vertexSubsetL = [1:nVertexL]';
    end
    
        
        nR = length(vertexSubsetR);
        nL = length(vertexSubsetL);
    
        subResidL = residL(:,vertexSubsetL);
        subResidR = residR(:,vertexSubsetR);
    
    %% gcv on left and right hemispheres with sandwich estimator:
    % this is done separately to obtain bandwidths used in whole-brain kernel
    % matrix
        

    outfull_swL = fullcovacem_sandwich(subResidL,sigmaemresults{k}.estsmsigmasqemL.sigmasqem,familyst,randLatL,randLongL,myhvec,1e-3,false,true);
    myresults.mseSA_L = outfull_swL.mseSA;
  
    myresults.hvecminL = outfull_swL.hvecmin;  
    myresults.mseSC_L = outfull_swL.mseSC;
    myresults.mseSEg_L = outfull_swL.mseSEg;
    
    fprintf('\nCompleted left hemisphere\n')
    
    % Edit 13 Febrauary: 
    %   bandwidth selection for left hemisphere only.
%     outfull_swR = fullcovacem_sandwich(residR,sigmaemresults.estsmsigmasqemR.sigmasqem,familyst,randLatR,randLongR,myhvec,1e-3,false,true);
%     myresults.hvecminR = outfull_swR.hvecmin;
%     myresults.mseSA_R = outfull_swR.mseSA;
%     myresults.mseSC_R = outfull_swR.mseSC;
%     myresults.mseSEg_R = outfull_swR.mseSEg;
%     myresults.hvec = outfull_swR.hvec;
%     clear outfull_swR;

    clear outfull_swL; 
    % calculate Sandwich estimator
    [kernL,~,unkernL] = createkernmat(randLatL,randLongL,mean(myresults.hvecminL),true);
    
    %[kernR,~,unkernR] = createkernmat(randLatR,randLongR,mean(myresults.hvecminR),true);
    [kernR,~,unkernR] = createkernmat(randLatR,randLongR,mean(myresults.hvecminL),true);
    
    kernLR = [kernL,  zeros(nL,nR); zeros(nR,nL), kernR];
    unkernLR = [unkernL, zeros(nL,nR); zeros(nR,nL), unkernR];
    clear kernL kernR unkernL unkernR; 
    
    % Sandwich estimate for whole cortex:  
    fprintf('\nBeginning Combined Left and Right Hemispheres SW Estimation\n')

    %myresults.outfull_swLR = fullcovacem_sandwich_inputkernmat([subResidL subResidR],...
    %    sigmaemresults{k}.smsigmasqemLR,familyst,kernLR,1000,false,1e-03);
    
    myresults.outfull_swLR = fullcovacem_sandwich_inputkernmat([subResidL subResidR],...
        sigmaemresults{k}.smsigmasqemLR,familyst,kernLR,neigsSA,neigsSC,neigsSEg,false,1e-03);
    
    %N - n1 = 1094 - 135
    
    
   fprintf('\nCompleted Combined Left and Right Hemispheres SW Estimation\n')

    % EIGS is much slower for a smallish number of vertices, e.g.: 2969: 
    % full eig is 9.64 sec; eigs with 1200 is 570 seconds!
    % 12 February 2017: changed from 1200 to 1000 due to convergence errors
    swresults{k} = myresults;
    
end

mytime = toc;
save('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_bandwidthselection_sw.mat',...
    'swresults','cvpartR','cvpartL','mytime','-v7.3');
