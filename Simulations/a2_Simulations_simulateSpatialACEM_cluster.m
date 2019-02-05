%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate the Spatial ACEM model --
%       ACE model with estimation of measurement error
%
%--------------------------------------------------------
 
tic;
% These simulations examine different estimates of the covariance
% matrix 
addpath './Functions'

% covariance matrices created in SimulateSpatialACE.m:
load './SupportingDataFiles/CovarianceMatricesForSimulationsSharmACEM_compress.mat'

% expand covariance matrices
sigmaa = sigmaa_v*sigmaa_l*sigmaa_v';
sigmac = sigmac_v*sigmac_l*sigmac_v';
sigmaeg = sigmaeg_v*sigmaeg_l*sigmaeg_v';


%% Key to indices:
%depthLabels = {'Sandwich Symm','Sandwich Trunc','SL Symm','SL Trunc','PSD Est Rank','PSD True Rank','MLE','SMLE','MWLE'};
depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};
nVertex = length(Lat);

outerloopniter = 21;
innerloopniter = 48; %e.g., equal to the number of workers.
% 21*48 = 1008 simulations in total

tniter = outerloopniter*innerloopniter;
simresultsCovfunMSE = zeros(3,6,outerloopniter*innerloopniter);
simresultsVarMSE = zeros(3,9,outerloopniter*innerloopniter);
simresultsCovfunAVE = zeros(3,9,nVertex,nVertex);
simresultsh2AVE = zeros(9,nVertex);
simresultsh2MSE = zeros(9,outerloopniter*innerloopniter);

% variance from both unique environmental variance and measurement error:
% calculated for MLE, SMMLE, MWLE:
simresultsSigmasqeAVE  = zeros(3,nVertex);
simresultsSigmasqeMSE = zeros(3,nVertex,outerloopniter*innerloopniter);

rng(714,'twister')
rngvector = randsample(10000,tniter);
allconvergence_estrank = zeros(tniter,1);
allconvergence_truerank = zeros(tniter,1);
gradnorm_estrank = zeros(tniter,1);
gradnorm_truerank = zeros(tniter,1);

gditer = 5000;

% EDIT THIS TO CUSTOMIZE FOR YOUR CLUSTER------>
% local_Copy.settings specifies the cluster set-up to allow the use of all
% logical processors (on my machine, 2x the number of physical processors)
% whereas matlab will usually only grab the physical processors
parallel.importProfile('~/Dropbox/SpatialACDE/Programs/local_Copy.settings')
poolobj = parpool('local_Copy',48);
poolobj.NumWorkers

% NOTE: There is an "outerloop" and "innerloop"
% where the innerloop stores all full covariances matrix estimates
% and hence requires the storage of 48 results when using 48 workers, 
% whereas the outerloop updates average covariance estimates for bias 
% calculations; this avoids storing all 1008 simulation results in memory,
% (instead storing 48 in inner loop + 1 that is the average)

for j=1:outerloopniter
    tindices = innerloopniter*(j-1)+1;
    tindicesm1 = tindices-1; %To meet parfor requirements
    simresultsCovfunAll = zeros(3,9,nVertex,nVertex,innerloopniter);
    simresultsh2All = zeros(9,nVertex,innerloopniter);
    simresultsSigmasqetAll = zeros(3,nVertex,innerloopniter);
    
    parfor t = tindices:tindices+innerloopniter-1
        rng(rngvector(t),'twister');

        %     %% Simulate data:
        simdata = simulatedataACEM(100,100,200,betas,sigmaa,sigmac,sigmaeg,sigmaem);

        %% Estimate Spatial ACE:
        % obtain estimates of unique environmental / measurement error:
        mleresults = allvertexmle_ACE(simdata.data, simdata.subjMat, simdata.familyst,false);

        smmle = gcv_smoothmle(mleresults,simdata.familyst,Lat,Long,[6:50],false,false);
        % Chooses moderate amount of smoothing for SigmaA and SigmaC, less
        % for SigmaE, this is consistent with the smoothness/roughness of
        % the respective covariance functions.
        
        % Use 5-fold cv in mwle as in SL FSEM.
        % calculates mse at a subset of locations to decrease computation
        % time; in these simulations, the smallest bandwidth is chosen
        select100 = zeros(nVertex,1);
        temp = randsample(nVertex,100);
        select100(temp) = 1;
        bw_cv5 = cv_ROI_mwle_ACE(mleresults,select100,[8:2:20],Lat,Long,simdata.familyst,false);  
        % chooses minimal smoothing
        
        mwleresults = allvertexmwle_ACE(mleresults,bw_cv5.hvecmin,Lat,Long,simdata.familyst,false);
        
        %5-fold cv implementation of measurement error estimate -- not used.
        %estsigmaem = estsigmasqem_cv(smmle.smresid,smmle.smsigmasqa,smmle.smsigmasqc,smmle.smsigmasqe,Lat,Long,[7.6:0.1:8.2],5);
    
        % GCV implementation:
        estsigmaem = estsigmasqem_gcv(smmle.smresid,simdata.familyst,smmle.smsigmasqa,...
            smmle.smsigmasqc,smmle.smsigmasqe,Lat,Long,[8:0.05:8.5]);
     
        outfull_sw = fullcovacem_sandwich(smmle.smresid,estsigmaem.sigmasqem,simdata.familyst,Lat,Long,[7.5:0.05:8.5],1e-3);
     
        %3. SL's estimator using bw from above:
        outsl = fullcovacem_sl_symm(smmle.smresid,simdata.familyst,Lat,Long,outfull_sw.hvecmin);

        %4. Constrained optimization:

        % automated rank selection of genetic covariance matrix:
        %mydA = estrank(outfull_sw.valSA,0.90);
        %mydA = min(sum(three<0.9),nSubject);
        % Fixed based on looking at multiple simulations: 
        mydA = 8;

        %automated rank selection of common environmental covariance matrix:
        %mydC = estrank(outfull_sw.valSC,0.90);
        %mydC = min(sum(three<0.9),nSubject);
        % Fixed based on looking at multiple simulations: 
        mydC = 8;

        % For SigmaSEg, clear break at 6:
        mydE = 6;
        
        outfull_psd_estrank = fullcovacem_con(smmle.smresid,estsigmaem.sigmasqem,mydA,mydC,mydE,simdata.familyst,...
            Lat,Long,outfull_sw.hvecmin,gditer,outfull_sw);
        allconvergence_estrank(t) = outfull_psd_estrank.convergence;
        gradnorm_estrank(t) = outfull_psd_estrank.gradnorm(end);

        outfull_psd_truerank = fullcovacem_con(smmle.smresid,smmle.smsigmasqe,5,5,6,simdata.familyst,Lat,Long,outfull_sw.hvecmin,gditer,outfull_sw);
        allconvergence_truerank(t) = outfull_psd_truerank.convergence;
        gradnorm_truerank(t) = outfull_psd_truerank.gradnorm(end);
        %-----------------------------------------------------------------
        covfunMSE = zeros(2,6);
        covfunMSE(1,1) = norm(sigmaa-outfull_sw.smSA_symm,'fro')^2;
        covfunMSE(1,2) = norm(sigmaa-outfull_sw.smSA_psd,'fro')^2;

        covfunMSE(1,3) = norm(sigmaa-outsl.smSA_symm,'fro')^2;
        covfunMSE(1,4) = norm(sigmaa-outsl.smSA_psd,'fro')^2;

        covfunMSE(1,5) = norm(sigmaa-outfull_psd_estrank.smSA_psd,'fro')^2;
        covfunMSE(1,6) = norm(sigmaa-outfull_psd_truerank.smSA_psd,'fro')^2;

        %-------------------------------------------
        covfunMSE(2,1) = norm(sigmac-outfull_sw.smSC_symm,'fro')^2;
        covfunMSE(2,2) = norm(sigmac-outfull_sw.smSC_psd,'fro')^2;

        covfunMSE(2,3) = norm(sigmac-outsl.smSC_symm,'fro')^2;
        covfunMSE(2,4) = norm(sigmac-outsl.smSC_psd,'fro')^2;

        covfunMSE(2,5) = norm(sigmac-outfull_psd_estrank.smSC_psd,'fro')^2;
        covfunMSE(2,6) = norm(sigmac-outfull_psd_truerank.smSC_psd,'fro')^2;

        %-------------------------------------------
        covfunMSE(3,1) = norm(sigmaeg-outfull_sw.smSEg_symm,'fro')^2;
        covfunMSE(3,2) = norm(sigmaeg-outfull_sw.smSEg_psd,'fro')^2;

        covfunMSE(3,3) = norm(sigmaeg-outsl.smSEg_symm,'fro')^2;
        covfunMSE(3,4) = norm(sigmaeg-outsl.smSEg_psd,'fro')^2;

        covfunMSE(3,5) = norm(sigmaeg-outfull_psd_estrank.smSEg_psd,'fro')^2;
        covfunMSE(3,6) = norm(sigmaeg-outfull_psd_truerank.smSEg_psd,'fro')^2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        varMSE = zeros(2,9);
        varMSE(1,1) = norm(diag(sigmaa-outfull_sw.smSA_symm),'fro')^2;
        varMSE(1,2) = norm(diag(sigmaa-outfull_sw.smSA_psd),'fro')^2;

        varMSE(1,3) = norm(diag(sigmaa-outsl.smSA_symm),'fro')^2;
        varMSE(1,4) = norm(diag(sigmaa-outsl.smSA_psd),'fro')^2;

        varMSE(1,5) = norm(diag(sigmaa-outfull_psd_estrank.smSA_psd),'fro')^2;
        varMSE(1,6) = norm(diag(sigmaa-outfull_psd_truerank.smSA_psd),'fro')^2;

        %-------------------------------------------
        varMSE(2,1) = norm(diag(sigmac-outfull_sw.smSC_symm),'fro')^2;
        varMSE(2,2) = norm(diag(sigmac-outfull_sw.smSC_psd),'fro')^2;

        varMSE(2,3) = norm(diag(sigmac-outsl.smSC_symm),'fro')^2;
        varMSE(2,4) = norm(diag(sigmac-outsl.smSC_psd),'fro')^2;

        varMSE(2,5) = norm(diag(sigmac-outfull_psd_estrank.smSC_psd),'fro')^2;
        varMSE(2,6) = norm(diag(sigmac-outfull_psd_truerank.smSC_psd),'fro')^2;

        varMSE(3,1) = norm(diag(sigmaeg-outfull_sw.smSEg_symm),'fro')^2;
        varMSE(3,2) = norm(diag(sigmaeg-outfull_sw.smSEg_psd),'fro')^2;

        varMSE(3,3) = norm(diag(sigmaeg-outsl.smSEg_symm),'fro')^2;
        varMSE(3,4) = norm(diag(sigmaeg-outsl.smSEg_psd),'fro')^2;

        varMSE(3,5) = norm(diag(sigmaeg-outfull_psd_estrank.smSEg_psd),'fro')^2;
        varMSE(3,6) = norm(diag(sigmaeg-outfull_psd_truerank.smSEg_psd),'fro')^2;

        % add MLE, SMLE and WMLE:
        varMSE(1,7) = sum((diag(sigmaa) - mleresults.sigmasqA').^2);
        varMSE(1,8) = sum((diag(sigmaa) - smmle.smsigmasqa').^2);
        varMSE(1,9) = sum((diag(sigmaa) - mwleresults.sigmasqA').^2);

        varMSE(2,7) = sum((diag(sigmac) - mleresults.sigmasqC').^2);
        varMSE(2,8) = sum((diag(sigmac) - smmle.smsigmasqc').^2);
        varMSE(2,9) = sum((diag(sigmac) - mwleresults.sigmasqC').^2);

        % unique genetic effects not identifiable in mle approaches
        % but calculated here for completeness:
        varMSE(3,7) = sum((diag(sigmaeg) - mleresults.sigmasqE').^2);
        varMSE(3,8) = sum((diag(sigmaeg) - smmle.smsigmasqe').^2);
        varMSE(3,9) = sum((diag(sigmaeg) - mwleresults.sigmasqE').^2);
        
        
        % Calculate ISE for h2:
        h2MSE = zeros(1,9);
        h2MSE(1) = sum((h2-outfull_sw.h2_symm).^2);
        h2MSE(2) = sum((h2-outfull_sw.h2_psd).^2);
        h2MSE(3) = sum((h2-outsl.h2_symm).^2);
        h2MSE(4) = sum((h2-outsl.h2_psd).^2);
        h2MSE(5) = sum((h2-outfull_psd_estrank.h2).^2);
        h2MSE(6) = sum((h2-outfull_psd_truerank.h2).^2);
        h2MSE(7) = sum((h2 - mleresults.h2).^2);
        h2MSE(8) = sum((h2 - smmle.h2).^2);
        h2MSE(9) = sum((h2 - mwleresults.h2).^2);

        % Calculate ISE for SigmasqEt:
        sigmasqetMSE = zeros(3,1);
        sigmaet = diag(sigmaeg)+sigmaem;
        sigmasqetMSE(1) = sum((sigmaet - mleresults.sigmasqE').^2);
        sigmasqetMSE(2) = sum((sigmaet - smmle.smsigmasqe').^2);
        sigmasqetMSE(3) = sum((sigmaet - mwleresults.sigmasqE').^2);

        simresultsCovfunMSE(:,:,t) = covfunMSE;
        simresultsVarMSE(:,:,t) = varMSE;
        simresultsh2MSE(:,t) = h2MSE;
        simresultsSigmasqeMSE(:,t) = sigmasqetMSE;
        
        % formatting for parfor:
        temp = zeros(3,6,nVertex,nVertex);
        temp(1,1,:,:) = reshape(outfull_sw.smSA_symm,[1,1,nVertex,nVertex]);
        temp(1,2,:,:) = reshape(outfull_sw.smSA_psd,[1,1,nVertex,nVertex]);
        temp(1,3,:,:) = reshape(outsl.smSA_symm,[1,1,nVertex,nVertex]);
        temp(1,4,:,:) = reshape(outsl.smSA_psd,[1,1,nVertex,nVertex]);
        temp(1,5,:,:) = reshape(outfull_psd_estrank.smSA_psd,[1,1,nVertex,nVertex]);
        temp(1,6,:,:) = reshape(outfull_psd_truerank.smSA_psd,[1,1,nVertex,nVertex]);
        temp(1,7,:,:) = reshape(diag(mleresults.sigmasqA),[1,1,nVertex,nVertex]);
        temp(1,8,:,:) = reshape(diag(smmle.smsigmasqa),[1,1,nVertex,nVertex]);
        temp(1,9,:,:) = reshape(diag(mwleresults.sigmasqA),[1,1,nVertex,nVertex]);

        temp(2,1,:,:) = reshape(outfull_sw.smSC_symm,[1,1,nVertex,nVertex]);
        temp(2,2,:,:) = reshape(outfull_sw.smSC_psd,[1,1,nVertex,nVertex]);
        temp(2,3,:,:) = reshape(outsl.smSC_symm,[1,1,nVertex,nVertex]);
        temp(2,4,:,:) = reshape(outsl.smSC_psd,[1,1,nVertex,nVertex]);
        temp(2,5,:,:) = reshape(outfull_psd_estrank.smSC_psd,[1,1,nVertex,nVertex]); 
        temp(2,6,:,:) = reshape(outfull_psd_truerank.smSC_psd,[1,1,nVertex,nVertex]); 
        temp(2,7,:,:) = reshape(diag(mleresults.sigmasqC),[1,1,nVertex,nVertex]);
        temp(2,8,:,:) = reshape(diag(smmle.smsigmasqc),[1,1,nVertex,nVertex]);
        temp(2,9,:,:) = reshape(diag(mwleresults.sigmasqC),[1,1,nVertex,nVertex]);

        temp(3,1,:,:) = reshape(outfull_sw.smSEg_symm,[1,1,nVertex,nVertex]);
        temp(3,2,:,:) = reshape(outfull_sw.smSEg_psd,[1,1,nVertex,nVertex]);
        temp(3,3,:,:) = reshape(outsl.smSEg_symm,[1,1,nVertex,nVertex]);
        temp(3,4,:,:) = reshape(outsl.smSEg_psd,[1,1,nVertex,nVertex]);
        temp(3,5,:,:) = reshape(outfull_psd_estrank.smSEg_psd,[1,1,nVertex,nVertex]); 
        temp(3,6,:,:) = reshape(outfull_psd_truerank.smSEg_psd,[1,1,nVertex,nVertex]); 
        temp(3,7,:,:) = reshape(diag(mleresults.sigmasqE),[1,1,nVertex,nVertex]);
        temp(3,8,:,:) = reshape(diag(smmle.smsigmasqe),[1,1,nVertex,nVertex]);
        temp(3,9,:,:) = reshape(diag(mwleresults.sigmasqE),[1,1,nVertex,nVertex]);

        simresultsCovfunAll(:,:,:,:,t-tindicesm1) = reshape(temp,[3,9,nVertex,nVertex,1]);
        
        
        temp = zeros(9,nVertex);
        temp(1,:) = outfull_sw.h2_symm;
        temp(2,:) = outfull_sw.h2_psd;
        temp(3,:) = outsl.h2_symm;
        temp(4,:) = outsl.h2_psd;
        temp(5,:) = outfull_psd_estrank.h2;
        temp(6,:) = outfull_psd_truerank.h2;
        temp(7,:) = mleresults.h2;
        temp(8,:) = smmle.h2;
        temp(9,:) = mwleresults.h2;
        simresultsh2All(:,:,t-tindicesm1) = reshape(temp,[9,nVertex,1]);
        
        % collect sigmaet averages:
        temp = zeros(3,nVertex);
        temp(1,:) = mleresults.sigmasqE;
        temp(2,:) = smmle.smsigmasqe;
        temp(3,:) = mwleresults.sigmasqE;
        simresultsSigmasqetAll(:,:,t-tindicesm1) = temp;
    end

    % average for bias calculations:
    temp = squeeze(mean(simresultsCovfunAll,5));
    simresultsCovfunAVE = simresultsCovfunAVE + temp./outerloopniter;

    temp = squeeze(mean(simresultsh2All,3));
    simresultsh2AVE = simresultsh2AVE + temp./outerloopniter;
    
    temp = mean(simresultsSigmasqetAll,3);
    simresultsSigmasqeAVE = simresultsSigmasqeAVE + temp./outerloopniter;
end

simresultsCovMatBias = zeros(3,9,nVertex,nVertex);
simresultsCovMatBiassq_sumV = zeros(3,6);

simresultsVarVecBias = zeros(3,9,nVertex);
simresultsVarVecBiassq_sumV = zeros(3,9);

simresultsh2Bias = zeros(9,nVertex);
simresultsh2Biassq_sumV = zeros(9,1);
for n=1:9
    simresultsCovMatBias(1,n,:,:) = squeeze(simresultsCovfunAVE(1,n,:,:)) - sigmaa;
    simresultsCovMatBiassq_sumV(1,n) = sum(sum(simresultsCovMatBias(1,n,:,:).^2));
    simresultsVarVecBias(1,n,:) = diag(squeeze(simresultsCovfunAVE(1,n,:,:))) - diag(sigmaa);
    simresultsVarVecBiassq_sumV(1,n) = sum(simresultsVarVecBias(1,n,:).^2);
    
    simresultsCovMatBias(2,n,:,:) = squeeze(simresultsCovfunAVE(2,n,:,:)) - sigmac;
    simresultsCovMatBiassq_sumV(2,n) = sum(sum(simresultsCovMatBias(2,n,:,:).^2));
    simresultsVarVecBias(2,n,:) = diag(squeeze(simresultsCovfunAVE(2,n,:,:))) - diag(sigmac);
    simresultsVarVecBiassq_sumV(2,n) = sum(simresultsVarVecBias(2,n,:).^2);
    
    simresultsCovMatBias(3,n,:,:) = squeeze(simresultsCovfunAVE(3,n,:,:)) - sigmaeg;
    simresultsCovMatBiassq_sumV(3,n) = sum(sum(simresultsCovMatBias(3,n,:,:).^2));
    simresultsVarVecBias(3,n,:) = diag(squeeze(simresultsCovfunAVE(3,n,:,:))) - diag(sigmaeg);
    simresultsVarVecBiassq_sumV(3,n) = sum(simresultsVarVecBias(3,n,:).^2);
    
    simresultsh2Bias(n,:) = h2 - simresultsh2AVE(n,:)';
    simresultsh2Biassq_sumV(n) = sum(simresultsh2Bias(n,:).^2);
end


simresultsCovMatMSE_sumV = squeeze(mean(simresultsCovfunMSE,3));
simresultsCovMatVariance_sumV = simresultsCovMatMSE_sumV - simresultsCovMatBiassq_sumV(:,1:6);
simresultsVarVecMSE_sumV = squeeze(mean(simresultsVarMSE,3));
simresultsVarVecVariance_sumV = simresultsVarVecMSE_sumV - simresultsVarVecBiassq_sumV;
simresultsh2MSE_sumV = mean(simresultsh2MSE,2);
simresultsh2Variance_sumV = simresultsh2MSE_sumV - simresultsh2Biassq_sumV;

save(['./aSimulations_100_100_200_' date],'sigmaa','sigmac','sigmaeg','sigmaem','betas',...
    'simresultsCovfunMSE','simresultsCovfunAVE','simresultsCovMatBias',...
    'simresultsVarMSE','simresultsVarVecBias',...
    'simresultsSigmasqeAVE','simresultsSigmasqeMSE',...
    'simresultsh2AVE','simresultsh2MSE','simresultsh2Bias',...
    'simresultsCovMatBiassq_sumV','simresultsCovMatVariance_sumV','simresultsCovMatMSE_sumV',...
    'simresultsVarVecBiassq_sumV','simresultsVarVecVariance_sumV','simresultsVarVecMSE_sumV',...
    'simresultsh2Biassq_sumV','simresultsh2Variance_sumV','simresultsh2MSE_sumV',...
    'rngvector','allconvergence_estrank','allconvergence_truerank','gradnorm_estrank','gradnorm_truerank','depthLabels');

toc
delete(poolobj)
