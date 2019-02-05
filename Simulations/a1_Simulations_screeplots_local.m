%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benjamin Risk
% Simulate the Spatial ACEM model --
%       ACE model with estimation of measurement error
% Create screeplots of eigenvalues to determine
% ranks used in simulations
%--------------------------------------------------------
 
% This code is used to determine the number of components
% in the PSD-ACE in subsequent simulations. Note that
% the subsequent programs can be run independently
% using the ranks 8, 8, and 6. 

tic;


addpath ./Functions
load ./Simulations/SupportingDataFiles/CovarianceMatricesForSimulationsSharmACEM_compress.mat

% expand covariance matrices
sigmaa = sigmaa_v*sigmaa_l*sigmaa_v';
sigmac = sigmac_v*sigmac_l*sigmac_v';
sigmaeg = sigmaeg_v*sigmaeg_l*sigmaeg_v';


%% Key to indices:
nVertex = length(Lat);
tniter = 20;
% takes a few minutes to run


rng(111,'twister')
rngvector = randsample(10000,tniter);

parallel.importProfile('~/Dropbox/SpatialACDE/Programs/local_Copy.settings')
poolobj = parpool('local_Copy',20);
poolobj.NumWorkers

allValSA = zeros(nVertex,tniter);
allValSC = zeros(nVertex,tniter);
allValSEg = zeros(nVertex,tniter);
tic;
parfor t = 1:tniter
        rng(rngvector(t),'twister');

        %% Simulate data:
        simdata = simulatedataACEM(100,100,200,betas,sigmaa,sigmac,sigmaeg,sigmaem);

        %% Estimate Spatial ACE:
        % obtain estimates of unique environmental / measurement error:
        mleresults = allvertexmle_ACE(simdata.data, simdata.subjMat, simdata.familyst,false);

        smmle = gcv_smoothmle(mleresults,simdata.familyst,Lat,Long,[6:60],false,false);
        % note: generally chooses a lot of smoothing, choice varies from simulation to simulation for sigmaa and sigmac; 
        
       % Fix bw at 8. <ExamineEstimationSigmaemSigmaeg>
       estsigmaem = estsigmasqem_gcv(smmle.smresid,simdata.familyst,smmle.smsigmasqa,...
            smmle.smsigmasqc,smmle.smsigmasqe,Lat,Long,8);
     
        
        outfull_sw = fullcovacem_sandwich(smmle.smresid,estsigmaem.sigmasqem,simdata.familyst,Lat,Long,[6:2:20],1e-3);
        allValSA(:,t) = outfull_sw.valSA;
        allValSC(:,t) = outfull_sw.valSC;
        allValSEg(:,t) = outfull_sw.valSEg;        
end
toc

save('./ExamineScreeplotsSimulateACEM.mat','allValSA','allValSC','allValSEg');


subplot(3,1,1)
plot(allValSA(1:20,:))
meanValSA = mean(allValSA,2);
line([8,8],[0,20],'col','red')
% 8 used in simulations


subplot(3,1,2)
plot(allValSC(1:20,:))
meanValSC = mean(allValSC,2);
line([8,8],[0,20],'col','red')
% 8 used in simulations


subplot(3,1,3)
plot(allValSEg(1:20,:))
line([6,6],[0,25],'col','red')
% 6 used in simulations


