%----------------------
% Estimate point-wise weighted MLE for cortical thickness
% Re-ran 15 October 2016 with 556 subjects
% Re-ran 11 April 2017 with 1200-subject release, 
%   after subsetting to twins and singletons, 676 subjects
%------------------------------------

load ~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/mle_cortThick676_HIGH.mat
load ~/Dropbox/SpatialACDE/Data/supportingdatafiles/latlong.mat
load '~/Dropbox/SpatialACDE/Data/supportingdatafiles/matrix_GordonRev2_CII_CORTEX_RIGHT.mat'
load ~/Dropbox/SpatialACDE/Data/subjectData_cortThickForMLE.mat
addpath '~/Dropbox/SpatialACDE/Programs/Functions'


%% Calculate WTD mle for left hemisphere:
[nSubject,nVertex] = size(mleresultsL_HIGH.fixefresid);

%% Calculate Vertex-wise weighted MLE:
sigmasq_a_wtdThick = zeros(3,nVertex);
sigmasq_n_wtdThick = zeros(2,nVertex);
lrt_a_wtdThick = zeros(nVertex,1);
pvalue_a_Thick = zeros(nVertex,1);
eFlag_a_wtdThick = lrt_a_wtdThick;
eFlag_n_wtdThick = eFlag_a_wtdThick;

% include vertices within h arclengths of focal vertex in non-zero weights.
% For example, for the first vertex,
% sum(dist<1) 1 (itself)
% sum(dist<2) 11
% sum(dist<5) 61
hvec = linspace(1,10,5);

%Estimate optimal bandwidth for MWLE using 5-fold cv on family id:
%To decrease computational expense, evaluate at approximately
%10% of the vertices:

rng(123)
randROI = binornd(1,0.10,nVertex,1);

rng(123)
tic;
cvout = cv_ROI_mwle_ACE(mleresultsL_HIGH, randROI, hvec, latL,longL, familyst);
toc
% two minutes    
figure;
plot(cvout.hvec,cvout.mse)
    
% redo with a finer search:
%hvec2 = linspace(1,3.5,10);
hvec2= [1.18,1.2:0.1:2.6];

%use same seed to get the same cv partition:
rng(123)
tic;
cvout2 = cv_ROI_mwle_ACE(mleresultsL_HIGH, randROI, hvec2, latL, longL, familyst);
toc

allhvec = [cvout.hvec,cvout2.hvec];
[allhvec,orderhvec] = sort(allhvec);
allmse = [cvout.mse;cvout2.mse];
allmse = allmse(orderhvec);

plot(allhvec,allmse)

a=plot(allhvec,allmse,':+','LineWidth',2);
xlabel('bandwidth (degrees)')
ylabel('MSE')
title('5-fold CV for weighted likelihood: Cortical Thickness')
saveas(a,'~/Dropbox/SpatialACDE/Documents/Figures/Figures_1200SubjectRelease/CV_MWLE_Thick_RandROI_10percent.png');

% h_ob: 2.3

h_ob = allhvec(allmse==min(allmse))

%%%%%%%%%%% Perform same analysis to see how the estimate changes:

rng(321)
randROI = binornd(1,0.10,nVertex,1);

rng(321)
cvout = cv_ROI_mwle_ACE(mleresultsL_HIGH, randROI, hvec, latL,longL, familyst);

%use same seed to get the same cv partition:
rng(321)
tic;
cvout2 = cv_ROI_mwle_ACE(mleresultsL_HIGH, randROI, hvec2, latL, longL, familyst);
toc

allhvec = [cvout.hvec,cvout2.hvec];
[allhvec,orderhvec] = sort(allhvec);
allmse = [cvout.mse;cvout2.mse];
allmse = allmse(orderhvec);

plot(allhvec,allmse)

a=plot(allhvec,allmse,':+','LineWidth',2);
xlabel('bandwidth (degrees)')
ylabel('MSE')
title('5-fold CV for weighted likelihood: Cortical Thickness')
saveas(a,'~/Dropbox/SpatialACDE/Documents/Figures/Figures_1200SubjectRelease/CV_MWLE_Thick_RandROI_10percent_version2.png');


allhvec(allmse==min(allmse))
% second set of locations chose 1.9



%% Try on full data: 
hvec2= [1.2:0.1:2.6];
randROI = ones(nVertexL,1);
%use same seed to get the same cv partition:
rng(321)
tic;
cvout2 = cv_ROI_mwle_ACE(mleresultsL_HIGH, randROI, hvec2, latL, longL, familyst);
toc


a=plot(cvout2.hvec,cvout2.mse,':+','LineWidth',2);
xlabel('bandwidth (degrees)')
ylabel('MSE')
title('5-fold CV for weighted likelihood: Cortical Thickness')
saveas(a,'~/Dropbox/SpatialACDE/Documents/Figures/Figures_1200SubjectRelease/CV_MWLE_Thick_RandROI_LeftCortex_AllVertices.png');

[b,c] = min(cvout2.mse)
hvec2(c)

h_ob = hvec2(c);


%% Inspect the selected bandwidth:
% on average, how many vertices are in the neighborhood for this optimal
% bandwidth:
randROI = binornd(1,0.10,nVertex,1);
indROI = find(randROI);
nneighbors = zeros(length(indROI),1);
parfor k = 1:sum(randROI)
     tindex = indROI(k);
     dist = double(distance(latL(tindex),longL(tindex),latL,longL));
     lindices = (dist<=h_ob);
     nneighbors(k) = sum(lindices);
end
var(nneighbors)
mean(nneighbors)
min(nneighbors)
max(nneighbors)


tic;
mwleresultsL_HIGH = allvertexmwle_ACE(mleresultsL_HIGH,h_ob,latL,longL,familyst);
toc

tabulate(mwleresultsL_HIGH.eFlag_a)
% <= 0 means it did not converge. 

scatter(mwleresultsL_HIGH.sigmasqA,mleresultsL_HIGH.sigmasqA)

tic;
mwleresultsR_HIGH = allvertexmwle_ACE(mleresultsR_HIGH,h_ob,latR,longR,familyst);
toc
tabulate(mwleresultsR_HIGH.eFlag_a)

save('~/Dropbox/SpatialACDE/Results/mwle_CortThick676_HIGH.mat','mwleresultsL_HIGH','mwleresultsR_HIGH')
