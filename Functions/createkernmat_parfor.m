function [ kernmat, rowsum, unkernmat ] = createkernmat_parfor(lat,long,h)
% Benjamin Risk
% 6 March 2017
%createkernmat Create smoothing matrix for local constant smoothing
%Input:
% lat
% long
% h - size of bandwidth

%Output:
% kernmat: V x V smoothing matrix normalized such that K1 = 1
% rowsum: the sum across columns of the kernel matrix before normalization
% unkernmat: not normalized version of kernmat

nVertex = length(lat);
unkernmat = zeros(nVertex,nVertex);


parfor v=1:nVertex
%   dists = distance(lat(v),long(v),lat,long);
    % 7 January 2016: Changed to use mygreatcirc since the cluster uses
    %matlab 2013a and negative values produce errors in distance()
    dists = mygreatcirc(lat(v),long(v),lat,long);
    temp = biweight(dists/h)/h;
    unkernmat(v,:) = temp;
%   kernmat(v,:) = biweight(dists/h)/h;
end
rowsum = sum(unkernmat,2);
kernmat = unkernmat.*repmat(1./rowsum,1,nVertex);
end

