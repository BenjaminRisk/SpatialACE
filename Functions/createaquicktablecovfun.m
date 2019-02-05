function [outtable] = createaquicktablecovfun(sigmaa,sigmac,sigmaeg,outsl,outfull_sw,outfull_psd_estrank,outfull_psd_truerank)
%CREATEAQUICKTABLECOVFUN Create a table that sums the squared differences between true
%covariance functions and estimates
%   Calculates errors between the truth and estimate from a single simulation, which is used for
%   illustration purposes. 
%   Calculates errors between the true heritability and estimate from a single simulation, which is used for
%   illustration purposes. 
% INPUT:
%   sigmaa -- true additive genetic covariance matrix, nvertices x
%       nvertices
%   sigmac -- true common environmental
%   sigmaeg -- true unique environmental
%   outsl -- estimates using FSEM (sl Shikai Luo)
%   outfull_sw -- estimates using sandwich variant of FSEM
%   outfull_psd_truerank -- PSD ACE estimate using true ranks
%   outfull_psd_estrank -- PSD ACE estimate using estimated ranks


        depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};
        covfunISE = zeros(2,6);
        covfunISE(1,1) = norm(sigmaa-outfull_sw.smSA_symm,'fro')^2;
        covfunISE(1,2) = norm(sigmaa-outfull_sw.smSA_psd,'fro')^2;

        covfunISE(1,3) = norm(sigmaa-outsl.smSA_symm,'fro')^2;
        covfunISE(1,4) = norm(sigmaa-outsl.smSA_psd,'fro')^2;

        covfunISE(1,5) = norm(sigmaa-outfull_psd_estrank.smSA_psd,'fro')^2;
        covfunISE(1,6) = norm(sigmaa-outfull_psd_truerank.smSA_psd,'fro')^2;

        %-------------------------------------------
        covfunISE(2,1) = norm(sigmac-outfull_sw.smSC_symm,'fro')^2;
        covfunISE(2,2) = norm(sigmac-outfull_sw.smSC_psd,'fro')^2;

        covfunISE(2,3) = norm(sigmac-outsl.smSC_symm,'fro')^2;
        covfunISE(2,4) = norm(sigmac-outsl.smSC_psd,'fro')^2;

        covfunISE(2,5) = norm(sigmac-outfull_psd_estrank.smSC_psd,'fro')^2;
        covfunISE(2,6) = norm(sigmac-outfull_psd_truerank.smSC_psd,'fro')^2;

        %-------------------------------------------
        covfunISE(3,1) = norm(sigmaeg-outfull_sw.smSEg_symm,'fro')^2;
        covfunISE(3,2) = norm(sigmaeg-outfull_sw.smSEg_psd,'fro')^2;

        covfunISE(3,3) = norm(sigmaeg-outsl.smSEg_symm,'fro')^2;
        covfunISE(3,4) = norm(sigmaeg-outsl.smSEg_psd,'fro')^2;

        covfunISE(3,5) = norm(sigmaeg-outfull_psd_estrank.smSEg_psd,'fro')^2;
        covfunISE(3,6) = norm(sigmaeg-outfull_psd_truerank.smSEg_psd,'fro')^2;
 
reSortIndexSix = [3,4,1,2,6,5];
SigmaA = covfunISE(1,reSortIndexSix)';
SigmaC = covfunISE(2,reSortIndexSix)';
SigmaEg = covfunISE(3,reSortIndexSix)';

outtable = table(SigmaA,SigmaC,SigmaEg,'RowNames',depthLabels(reSortIndexSix));

end

