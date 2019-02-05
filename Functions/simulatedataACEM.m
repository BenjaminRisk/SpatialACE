function [sim] = simulatedataACEM(n1,n2,n3,betas,sigmaa,sigmac,sigmaeg,sigmaem)
%simulatedataACEM
% Simulate Spatial ACE model with measurement error
%
%INPUT:
% n1: number of monozygotic twin pairs
% n2: number of dizygotic twin pairs
% n3: number of singletons
% sigmaa: V x V covariance matrix of additive genetic effects, where V is the number of locations
% sigmac: V x V covariance matrix of common environmental effects, where V
%           is the number of locations
% sigmaeg: V x V covariance matrix of unique environmental variance (also a Gaussian process)  
% sigmael: V x 1 measurement error (independent across space) 
%
%OUTPUT:
% sim is a structure with the following fields
% .familyst  Family structure. 
%       .MZtp1    indicator variable for individuals that are monozygotic twin pair one
%       .MZtp2    indicator variable for individuals that are monozygotic twin
%                   pair two
%       .DZtp1    indicator variable for dizygotic twin pair one
%       .DZtp2
%       .MDti     indicator variable for singletons. 
%       .familyID numeric identifier of family
%
% .aIJ  additive genetic effects of dizygotic twins that are unique to an
%           individual; see manuscript for an explanation.
% .aI   additive genetic effects shared by twins.
% .cI   common environmental effects
% .egIJ unique environmental effects
%   Note: measurement error is included in ``data''; see line 67.
% .data response matrix, N x V where N is the number of subjects
% .subjMat covariates;  N x 2 matrix where first column is a column of ones
%   and the second column is simulated from standard normal
%% Define groups and create design matrix
nind1 = n1*2;
nind2 = n2*2;
nVertex = size(sigmac,1);
N = nind1 + nind2 + n3;
DZ = zeros(N,1);
DZ(nind1+1:nind2+nind1) = 1;
subjMat = [ones(nind1+nind2+n3,1), randn(nind1+nind2+n3,1)];
sim.familyst.MZtp1 = logical([repmat([1;0],n1,1);zeros(nind2,1);zeros(n3,1)]);
sim.familyst.MZtp2 = logical([repmat([0;1],n1,1);zeros(nind2,1);zeros(n3,1)]);
sim.familyst.DZtp1 = logical([zeros(nind1,1);repmat([1;0],n2,1);zeros(n3,1)]);
sim.familyst.DZtp2 = logical([zeros(nind1,1);repmat([0;1],n2,1);zeros(n3,1)]);
sim.familyst.MDti = logical([zeros(nind1+nind2,1);ones(n3,1)]);

% create family id:
sim.familyst.familyID = [kron([1:(n1+n2)],[1 1]),(n1+n2+1):(n1+n2+n3)];


%% Simulate Spatial ACE:
rankA = rank(sigmaa);
rankC = rank(sigmac);
rankE = rank(sigmaeg);
sim.aIJ = [zeros(nVertex,nind1),simGP(sigmaa,rankA,nind2),zeros(nVertex,n3)]';
sim.aI = kron(simGP(sigmaa,rankA,n1+n2),[1 1]);
sim.aI = [sim.aI,simGP(sigmaa,rankA,n3)]';
sim.cI = kron(simGP(sigmac,rankC,n1+n2),[1 1]);
sim.cI = [sim.cI,simGP(sigmac,rankC,n3)]';
sim.egIJ = simGP(sigmaeg,rankE,N)';

sim.data = zeros(N,nVertex);
for v=1:nVertex
    sim.data(:,v) = subjMat*betas(v,:)' + sqrt(0.5)*(DZ.*(sim.aIJ(:,v)+sim.aI(:,v))) + (1-DZ).*sim.aI(:,v)...
        + sim.cI(:,v) + sim.egIJ(:,v) + sqrt(sigmaem(v))*randn(N,1);
end
sim.subjMat = subjMat;
end
