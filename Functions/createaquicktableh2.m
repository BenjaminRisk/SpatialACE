function [ outtable ] = createaquicktableh2(h2,outsl,outfull_sw,outfull_psd_truerank,outfull_psd_estrank,mleresults,smmle,mwleresults)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
     % Calculate MSE for h2:
        depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};
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

  

reSortIndexNine = [7,9,8,3,4,1,2,6,5];
h2MSE = h2MSE(reSortIndexNine)';
outtable = table(h2MSE,'RowNames',depthLabels(reSortIndexNine));      

end

