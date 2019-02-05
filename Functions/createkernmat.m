function [ kernmat, rowsum, unkernmat ] = createkernmat(lat,long,h,usesparse)
%createkernmat Create smoothing matrix for local constant smoothing
%Input:
% lat
% long
% h - size of bandwidth
% usesparse - true or false indicating whether to use sparse matrices
%Output:
% kernmat: V x V smoothing matrix normalized such that K1 = 1
% rowsum: the sum across columns of the kernel matrix before normalization
% unkernmat: not normalized version of kernmat
nVertex = length(lat);
if nargin<4
    usesparse = true;
end
if usesparse
    unkernmat = sparse(nVertex,nVertex);
else
    unkernmat = zeros(nVertex,nVertex);
end

for v=1:nVertex
    %dists = distance(lat(v),long(v),lat(v:end),long(v:end));
    % 7 January 2016: Changed to use mygreatcirc since the cluster uses
    %matlab 2013a and negative values produce errors in distance()
    dists = mygreatcirc(lat(v),long(v),lat(v:end),long(v:end));
    indexs = find(dists<=h);
    temp = biweight(dists(indexs)/h)/h;
    unkernmat(v,indexs+v-1) = temp;
    unkernmat(indexs+v-1,v) = temp;
%   kernmat(v,:) = biweight(dists/h)/h;
end
% using indices is slightly faster
% 1 minute to calculate the kernel
rowsum = sum(unkernmat,2);
kernmat = unkernmat.*repmat(1./rowsum,1,nVertex);
end

