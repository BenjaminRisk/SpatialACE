function [mydata] = readincortexst(indir,hemisphere)
% This function extracts the specified hemisphere from the cifti file
% and converts to a matlab matrix. The data from session 1 (RL phase 
% encoding) and session 2 (LR phase encoding) are concatenated.

hemisphere = validatestring(hemisphere, {'CORTEX_RIGHT','CORTEX_LEFT'});
side = 'R';
if strcmp(hemisphere,'CORTEX_LEFT')
    side = 'L';
end

% NOTE, if wb_command produces an error, try editing the lines below 
% to contain the full path to wb_command. 
unix(['wb_command -cifti-separate ' indirRL ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_RL.cortex.' side '.func.gii']);

% Example where the full path to wb_command is specified:
%unix(['/usr/bin/wb_command -cifti-separate ' indir ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_cortex.' side '.shape.gii']);

giidata = gifti([tempdir 'temp_cortex.' side '.shape.gii']);

mydata = giidata.cdata;

unix(['rm ' tempdir  'temp_cortex.' side '.shape.gii']);

load(['~/Dropbox/SpatialACDE/Data/supportingdatafiles/mapGIFTItoCIFTI_cortex_',side,'.mat']);
if strcmp(hemisphere,'CORTEX_RIGHT')
    mydata = double(mydata(mapping_CORTEX_RIGHT,:));
else mydata = double(mydata(mapping_CORTEX_LEFT,:));
end
