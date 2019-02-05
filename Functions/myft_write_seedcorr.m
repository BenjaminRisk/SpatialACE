function [matout] = myft_write_seedcorr(seedlist,psdresults,outfile,pathtotemplate,pathtomapGIFTItoCIFTI_L,pathtomapGIFTItoCIFTI_R)

% Calculate the covariance function of genetic component for the locations
% indexed by seedlist
% Input: 
%   seedlist:  (SURFACE VERTEX index + 1) for CortexLeft, i.e., where
%       SURFACE VERTEX is Gifti index using "Identify Surface Vertex" and
%       CortexLeft in wb_view
%       (SURFACE VERTEX index + 32492) for CortexRight
%   psdresults: a cell containing the divide and conquer solutions with
%       objects psdresults{k}.smXa (59412 x p_A) and psdresults{k}.smXc (59412
%       x p_C)
%   outfile: name of output; see defaults
%   pathtotemplate: path to any .dtseries.nii file in cifti format, used as a template
%   pathtomapGIFTItoCIFTI_L: path to the file /mapGIFTItoCIFTI_cortex_L
%   pathtomapGIFTItoCIFTI_R: path to the file /mapGIFTItoCIFTI_cortex_R
% matout: matrix with cifti dimensions. typically not output. 

if nargin<3
    outfile = ['~/Dropbox/SpatialACDE/Results/CorrACEMSigmaA_conpsd_seeds',num2str(min(seedlist)),'to',num2str(max(seedlist))];
end

if nargin<4
	test = ft_read_cifti('~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii');
else
	test = ft_read_cifti(pathtotemplate);
end

if nargin<5 
	load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L
	load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R
else
	load(pathtomapGIFTItoCIFTI_L)
	load(pathtomapGIFTItoCIFTI_R)
end

K = size(psdresults,1);

nt = length(seedlist);
test.time = seedlist;
%test.time = 1:nt;

test.dtseries = zeros(96854,nt);
rankA = size(psdresults{1}.smXa,2);

for k=1:K
    input1 = psdresults{k}.smXa;
    nanindices = isnan(input1(:,1));
    if sum(nanindices)
        if k==1
            input1(nanindices,:) = psdresults{2}.smXa(nanindices,:);
        elseif k==2 
            input1(nanindices,:) = psdresults{1}.smXa(nanindices,:);
        else 
            warning('There are missing data. Code must be modified for K>2.')
        end
    end
    gii_L = zeros(32492,rankA);
    gii_L(mapping_CORTEX_LEFT,:) = input1(1:length(mapping_CORTEX_LEFT),:);
    gii_R = zeros(32492,rankA);
    gii_R(mapping_CORTEX_RIGHT,:) = input1(length(mapping_CORTEX_LEFT)+1:end,:);
    gii_LR = [gii_L;gii_R];
    input2 = gii_LR(seedlist,:);
    input = gii_LR*input2';
    varAll = sum(gii_LR.^2,2);
    varSeed = sum(input2.^2,2);
    temp1 = sparse(1:64984,1:64984,1./sqrt(varAll));
    temp2 = diag(1./sqrt(varSeed));
    input = temp1*input*temp2;
    test.dtseries(1:64984,:) = test.dtseries(1:64984,:)+input./K;
end

ft_write_cifti(outfile,test,'parameter','dtseries');
matout = test.dtseries;

end
