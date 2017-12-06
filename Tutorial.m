%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUTORIAL.M
% 6 March 2017
% Benjamin Risk
% 
% PURPOSE: 
%   Fit SpatialACE with measurement error to simulated data.
%   Data are simulated on a sphere and the kernels use
%   geodesic distance via the great circle formula.
%
% Questions or comments:
%     email Ben at benjamin.risk@emory.edu
%--------------------------------------------------------

% These simulations examine different estimates of the covariance
% matrix 

% REVISE to match your directory structure:
cd ~/Dropbox/SpatialACDE/Package


addpath './Functions'

% covariance matrices constructed from covariance functions with spherical harmonics:
% 	covariance matrices constructed in <SimulateSpatialACE.m>, which uses
%	other scripts from Matlab central and is not provide here. Contact
%	benjamin.risk@emory.edu for additional information.
load './SupportingDataFiles/CovarianceMatricesForSimulationsSharmACEM.mat'


%% Key to indices:
depthLabels = {'S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE'};
nVertex = length(Lat);

rng(111,'twister')

%% Simulate data:
% simulate data for 100 MZ pairs, 100 DZ pairs, 200 singletons
% NOTE: sigmaem corresponds to measurement error, i.e., \sigma_{e,L}^2 in
% manuscript. 
simdata = simulatedataACEM(100,100,200,betas,sigmaa,sigmac,sigmaeg,sigmaem);

%% Estimate Spatial ACE:

%% STEP 1:  obtain estimates of unique environmental / measurement error:
tic;
mleresults = allvertexmle_ACE(simdata.data, simdata.subjMat, simdata.familyst,true);
toc
% 17 seconds with parfor

%% STEP 2: smooth (regularize) variance and covariate estimates:
smmle = gcv_smoothmle(mleresults,simdata.familyst,Lat,Long,[6:50],false,true);

%% STEP 3: estimate measurement error
estsigmaem = estsigmasqem_gcv(smmle.smresid,simdata.familyst,smmle.smsigmasqa,...
            smmle.smsigmasqc,smmle.smsigmasqe,Lat,Long,[8:0.05:8.5]);

%% STEP 4: stimate bandwidths and initial values using sandwich estimator:
outfull_sw = fullcovacem_sandwich(smmle.smresid,estsigmaem.sigmasqem,simdata.familyst,Lat,Long,[7.5:0.05:8.5],1e-3);        

%% STEP 5: estimate rank of covariance functions:
subplot(2,2,1)
plot(outfull_sw.valSA(1:40))
line([8,8],[0,16],'col','red')

% plot(log(outfull_sw.valSA(1:40)))
% choose rank; err on the SMALLER rank, but results are robust either way.     
mydA = 8; % 8 is a number chosen based on averaging screeplots from hundreds of simulations;
          % you can vary this number and see the results do not change very
          % much. Larger values tend to add a little bias. 

%generate estimates of rank of common environmental covariance:
subplot(2,2,2)
plot(outfull_sw.valSC(1:40))
line([6,6],[0,15],'col','red')
line([8,8],[0,15],'col','green')
mydC = 8; % see note from above; in the manuscript, used 8

% For SigmaSEg, clear break at 6:
subplot(2,2,3)
plot(outfull_sw.valSEg(1:40))
line([6,6],[0,25],'col','red')
mydE = 6; %unique environmental clearly 6
        
    
%% STEP 6: Estimate covariance functions
gditer = 5000;
tic;
outfull_psd_estrank = fullcovacem_con(smmle.smresid,estsigmaem.sigmasqem,mydA,mydC,mydE,simdata.familyst,...
            Lat,Long,outfull_sw.hvecmin,gditer,outfull_sw);
toc
% Takes about 30 seconds to run.

outfull_psd_truerank = fullcovacem_con(smmle.smresid,smmle.smsigmasqe,5,5,6,simdata.familyst,Lat,Long,outfull_sw.hvecmin,gditer,outfull_sw);
% You may receive many messages about "Decreasing lambda." This sometimes occurs even though
% the algorithm has adequately converged. The parameters in fullcovacem_con
% are relatively conservative.
% If you receive " Warning: size of gradient increased for vanishingly small lambda." Look at the message output 
% ``Size relative to initial gradient'', and you can judge whether
% the gradient is small enough. Size = ||current gradient|| / ||initial gradient||
% Increasing rank of bases can improve convergence. 



%% OPTIONAL: 
% Examine other methods
select100 = zeros(nVertex,1);temp = randsample(nVertex,100);select100(temp) = 1;

bw_cv5 = cv_ROI_mwle_ACE(mleresults,select100,[8:2:20],Lat,Long,simdata.familyst,false);  
mwleresults = allvertexmwle_ACE(mleresults,bw_cv5.hvecmin,Lat,Long,simdata.familyst,false);

%Luo et al. estimator using bw from above:
outsl = fullcovacem_sl_symm(smmle.smresid,simdata.familyst,Lat,Long,outfull_sw.hvecmin);



%% Example plots:
rng(321)
iIndex = randsample(1002,1);

% create a seed plot of additive genetic covariance from the vertex iIndex:
mypos =  [-22.9690   -0.2004   19.2740];
createaquickseedplotsigmaa(iIndex,sigmaa,outsl,outfull_sw,outfull_psd_estrank,outfull_psd_truerank,mypos,face,x,y,z,depthLabels);
  
% a different location
createaquickseedplotsigmaa(1,sigmaa,outsl,outfull_sw,outfull_psd_estrank,outfull_psd_truerank,mypos,face,x,y,z,depthLabels);
        
%-----------------------------------------------------------------
%% COMPARE MSEs:
% covariance functions:
createaquicktablecovfun(sigmaa,sigmac,sigmaeg,outsl,outfull_sw,outfull_psd_truerank,outfull_psd_estrank)

% heritability:
createaquicktableh2(h2,outsl,outfull_sw,outfull_psd_truerank,outfull_psd_estrank,mleresults,smmle,mwleresults)


      
