%----------------------
% Estimate cortical thickness heritability
% using the ACEM model with measurement error
% Updated 17 April 2017:
%   Re-ren for 1200-subject release
%------------------------------------

tic;
addpath ~/Dropbox/SpatialACDE/Programs/Functions
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/smmle_cortThick676_HIGH.mat
load ~/Dropbox/SpatialACDE/Data/subjectData_cortThickForCovfun.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/SigmaemCortThickACEM.mat
load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_bandwidthselection_sw.mat

K = 1; % number of data partitions
nsteps=1000; % number of steps in gradient descent algorithm
ntwins = sum(familyst.MZtp1)+sum(familyst.DZtp1);
nSubjects = size(dataMatL,2);
rankA = ntwins;
rankC = ntwins;
rankEg = nSubjects - sum(familyst.MZtp1);
nVertexR = size(dataMatR,1);
nVertexL = size(dataMatL,1);

% NOTE: THIS SECTION MODIFIED FOR DEBUGGING ONLY!!!
%for k=1:K
%     swresults{k}.hvecminL = 5;
%     swresults{k}.hvecminR = 5;
% end

% temporary code for NaNs that result when estimating too many eigenvalues:
for k=1:K
    yesA = ~isnan(swresults{k}.outfull_swLR.valSA);
    swresults{k}.outfull_swLR.valSA = swresults{k}.outfull_swLR.valSA(yesA);
    swresults{k}.outfull_swLR.vecSA = swresults{k}.outfull_swLR.vecSA(:,yesA);
    yesC = ~isnan(swresults{k}.outfull_swLR.valSC);
    swresults{k}.outfull_swLR.valSC = swresults{k}.outfull_swLR.valSC(yesC);
    swresults{k}.outfull_swLR.vecSC = swresults{k}.outfull_swLR.vecSC(:,yesC);
end

maxNumCompThreads(24);

psdresults = cell(K,1);
%% Loop through data partitions
for k=1:K
    if k>1
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
    
    residL = smmleL_HIGH.smresid(:,vertexSubsetL);    
    residR = smmleR_HIGH.smresid(:,vertexSubsetR);
    mybwL = mean(swresults{k}.hvecminL);
    %mybwR = mean(swresults{k}.hvecminR);
    

    %3. Constrained optimization:
    [kernL,~,unkernL] = createkernmat(randLatL,randLongL,mybwL,true);
    %[kernR,~,unkernR] = createkernmat(randLatR,randLongR,mybwR,true);
    [kernR,~,unkernR] = createkernmat(randLatR,randLongR,mybwL,true);
    
    unkernLR = [unkernL, zeros(nL,nR); zeros(nR,nL), unkernR];
    clear unkernL unkernR

    outfull_psdLR = fullcovacem_con_inputunkernmat([residL residR],...
            sigmaemresults{k}.smsigmasqemLR, rankA,rankC,rankEg,familyst,unkernLR,nsteps,swresults{k}.outfull_swLR,false);
        
  if K>1      
    %% CURRENT APPROACH: Use a robust inverse that removes small eigenvalues (which blow up in inverse)
    %kernL = createkernmat(randLatL,randLongL,mybwL,false);
    [ukernL,dkernL,vkernL] = svd(full(kernL)); 
    dkernL = diag(dkernL);

    % nsv = 5000;
    % max eigenvalue approximately equals 1
    % larger bw = faster decay
    mytolrobust = 0.0001;
    nsv = find(dkernL>mytolrobust,1,'last');
    robustInverseL = vkernL(:,1:nsv)*diag(1./dkernL(1:nsv))*ukernL(:,1:nsv)';

    %kernR = createkernmat(randLatR,randLongR,mybwR,false);
    [ukernR,dkernR,vkernR] = svd(full(kernR));
    dkernR = diag(dkernR);
    nsv = find(dkernR>mytolrobust,1,'last');
    robustInverseR = vkernR(:,1:nsv)*diag(1./dkernR(1:nsv))*ukernR(:,1:nsv)';

    myresults.robustunsmXa = [robustInverseL*outfull_psdLR.Xa(1:nL,:);robustInverseR*outfull_psdLR.Xa(nL+1:end,:)];
    myresults.robustunsmXc = [robustInverseL*outfull_psdLR.Xc(1:nL,:);robustInverseR*outfull_psdLR.Xc(nL+1:end,:)];
    myresults.robustunsmXeg = [robustInverseL*outfull_psdLR.Xeg(1:nL,:);robustInverseR*outfull_psdLR.Xeg(nL+1:end,:)];
        %     tempsmXa = [kernL, zeros(nL,nR);zeros(nR,nL), kernR]*myresults.robustunsmXa;
    %     norm(tempsmXa*tempsmXa' - outfull_psdLR.Xa*outfull_psdLR.Xa','fro')
    % 
    %     tempsmXc = [kernL, zeros(nSubvertexL,nSubvertexR);zeros(nVertexL,nVertexR), kernR]*robustunsmXc;
    %     norm(tempsmXc*tempsmXc' - outfull_psdLR.Xc*outfull_psdLR.Xc','fro')

    %     robustsmXa = rownormRbasisLR*robustunsmXa;
    %     robustsmXc = rownormRbasisLR*robustunsmXc;

    %     robustsampleFullLR_SA = robustsmXa*robustsmXa(startVertex:startVertex+1000,:)';

  end 
    myresults.outfull_psdLR = outfull_psdLR;
    psdresults{k} = myresults;   
end

%save('~/Dropbox/SpatialACDE/Results/AllEstimatesCortThick_subset_rkhs.mat','mleresultsL','mleresultsR',...
%    'smmleL','smmleR','outfull_psdLR','outfull_swLR','randL','randR','randLatL','randLatR',...
%    'surfdata_R','surfdata_L','mybwL','mybwR','lowDataL','lowDataR','unsmXa','unsmXc','-v7.3')
mytime = toc

save(['~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/AllEstimatesCortThickACEM_divideandconquer_K=' num2str(K) '.mat'],...
    'psdresults','mytime','-v7.3');


