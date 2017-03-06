function [outtable] = createaquicktablecovfun(sigmaa,sigmac,sigmaeg,outsl,outfull_sw,outfull_psd_estrank,outfull_psd_truerank)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
        depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};
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
 
reSortIndexSix = [3,4,1,2,6,5];
SigmaA = covfunMSE(1,reSortIndexSix)';
SigmaC = covfunMSE(2,reSortIndexSix)';
SigmaEg = covfunMSE(3,reSortIndexSix)';

outtable = table(SigmaA,SigmaC,SigmaEg,'RowNames',depthLabels(reSortIndexSix));

end

