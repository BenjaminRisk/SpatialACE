function [eigvec,eigval] = eig_descend(symmat,useeigs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[eigvec,eigval] = eig(symmat,'vector');
% eig has different behaviors in different versions of matlab!!!
% argument 'vector' not allowed on cluster in matlab 2013b
% useeigs: integer from 0 to V determining whether to use eigs
%   if no value is provided, then eig is used.
% 
% NOTE: if matrix is not symmetric, then the eigenvalues
% are complex and are sorted by their magnitudes!
% Returns in error if this occurs
if nargin<2
    useeigs=0;
end

if useeigs
    [eigvec,eigval] = eigs(symmat,useeigs,'la');
    eigval = diag(eigval);
else
    [eigvec,eigval] = eig(symmat);
    eigval = diag(eigval);
end

if ~isreal(eigval)
    error('Contains complex eigenvalues')
end

%re-order from largest to smallest:
[eigval,tindex] = sort(eigval,'descend');
eigvec = eigvec(:,tindex);

end

