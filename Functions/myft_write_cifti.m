function [] = myft_write_cifti(input,outfile,pathtotemplate,pathtomapGIFTItoCIFTI_L,pathtomapGIFTItoCIFTI_R)
% Create a .dtseries.nii file for viewing in wb_view
%INPUT:
% input: name of 59,412 x p matrix
% outfile: name with directory structure for output cifti file .dt.series.nii will be appended
% pathtotemplate: path to any .dtseries.nii file in cifti format, used as a template
% pathtomapGIFTItoCIFTI_L: path to the file /mapGIFTItoCIFTI_cortex_L
% pathtomapGIFTItoCIFTI_L: path to the file /mapGIFTItoCIFTI_cortex_R

%OUTPUT
% .dtseries.nii file for viewing in wb_view


if nargin<3
	test = ft_read_cifti('~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii');
else
	test = ft_read_cifti(pathtotemplate);
end

if nargin<4 
	load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_L
	load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_R
else
	load(pathtomapGIFTItoCIFTI_L)
	load(pathtomapGIFTItoCIFTI_R)
end

nt = size(input,2);

test.time = 1:nt;
test.dtseries = zeros(96854,nt);
gii_L = zeros(32492,nt);
gii_L(mapping_CORTEX_LEFT,:) = input(1:length(mapping_CORTEX_LEFT),:);
gii_R = zeros(32492,nt);
gii_R(mapping_CORTEX_RIGHT,:) = input(length(mapping_CORTEX_LEFT)+1:end,:);

test.dtseries(1:32492,:) = gii_L;
test.dtseries(32493:64984,:) = gii_R;

ft_write_cifti(outfile,test,'parameter','dtseries');
end
