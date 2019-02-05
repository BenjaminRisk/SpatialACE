function [svecSA,valSA,svecSC,valSC, svecSCnull, valSCnull] = fastcovace(R,MZtp1,MZtp2,DZtp1,DZtp2,MDti,lat,long,hvec,calcnull,usesparse)
%FASTCOVACE calculates a decomposition of the smoothed VxV covariance matrices
%   for the additive and common environmental variance components. 
%   This function uses random projections to avoid calculating the V x V
%   residual matrices and thus is scalable to large covariance matrices.
%
%INPUT:
% R: N x V matrix; residual matrix
% MZtp1
% MZtp2
% DZtp1
% DZtp2
% MDti
% lat V x 1 vector of latitudes
% long V x 1 vector of longitudes
% h: scalar; bandwidth
% calcnull: calculate SigmaC under the null hypothesis that SigmaA = 0
% usesparse: use sparse matrices; slower but scalable to massive matrices

%OUTPUT:
% svecSA: V x d matrix of eigenvectors for the additive genetic component, where d <= N
% svalSA: d x d matrix of associated positive eigenvalues
% svecSC: V x d matrix of eigenvectors for the common environmental
%           component
% svalSC: d x d matrix of associated eigenvalues
% svecSCnull: V x d matrix of eigenvectors under null hypothesis that SA=0
% svalSCnull: d x d matrix of associated eigenvalues


% create kernel matrix:
for k=1:length(bw)
    kernmat = createkernmat(lat,long,h,usesparse);
    [svecSA,valSA,svecSC,valSC, svecSCnull, valSCnull] = ifastcovace(R,MZtp1,MZtp2,DZtp1,DZtp2,MDti,kernmat,calcnull);
    %calculate elements of covariance matrix in blocks:
    
    
end
end
