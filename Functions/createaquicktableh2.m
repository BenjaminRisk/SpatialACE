function [ outtable ] = createaquicktableh2(h2,outsl,outfull_sw,outfull_psd_truerank,outfull_psd_estrank,mleresults,smmle,mwleresults)
%CREATE Summary of this function goes here
%   Calculates errors between the true heritability and estimate from a single simulation, which is used for
%   illustration purposes. 
% INPUT:
%   h2 -- true heritability
%   outsl -- estimates using FSEM (sl Shikai Luo)
%   outfull_sw -- estimates using sandwich variant of FSEM
%   outfull_psd_truerank -- PSD ACE estimate using true ranks
%   outfull_psd_estrank -- PSD ACE estimate using estimated ranks
%   mleresults -- MLE estimates
%   smmle -- smoothed MLE estimates
%   mwleresults -- maximum weighted likelihood as in SL FSEM


     % Calculate ISE for h2:
        depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};
        h2ISE = zeros(1,9);
        h2ISE(1) = sum((h2-outfull_sw.h2_symm).^2);
        h2ISE(2) = sum((h2-outfull_sw.h2_psd).^2);
        h2ISE(3) = sum((h2-outsl.h2_symm).^2);
        h2ISE(4) = sum((h2-outsl.h2_psd).^2);
        h2ISE(5) = sum((h2-outfull_psd_estrank.h2).^2);
        h2ISE(6) = sum((h2-outfull_psd_truerank.h2).^2);
        h2ISE(7) = sum((h2 - mleresults.h2).^2);
        h2ISE(8) = sum((h2 - smmle.h2).^2);
        h2ISE(9) = sum((h2 - mwleresults.h2).^2);

  

reSortIndexNine = [7,9,8,3,4,1,2,6,5];
h2ISE = h2ISE(reSortIndexNine)';
outtable = table(h2ISE,'RowNames',depthLabels(reSortIndexNine));      

end

