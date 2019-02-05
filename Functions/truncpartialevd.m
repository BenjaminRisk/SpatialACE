function [eigvec,eigval] = truncpartialevd(symmat,n)
% truncpartialevd Calculate the first n eigenvalues for a symmetric matrix,
% then truncate to positive eigenvalues
%INPUT:
% symmat    a symmetric matrix
% n         number of eigenvalue/eigenvectors pairs to extract
%OUTPUT: 
% eigvec    eigenvectors associated with positive eigenvalues
% eigval    diagonal matrix of positive eigenvalues


[vec,val] = eigs(symmat,n,'la');
val = diag(val);

% matlab notes that when using eigs, eigenvalues may not be sorted correctly:
[~,sindex] = sort(val,'descend');
vec = vec(:,sindex);
val = val(sindex);
tol = val(1)*0.01;
index = val>tol;
if sum(index)==length(sindex) 
    warning('Use a larger n -- all estimated eigenvalues are greater than 0.01*(largest eigenvalue)')
end
     
eigvec = vec(:,index);
eigval = diag(val(index));

end
