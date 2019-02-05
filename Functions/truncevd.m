function [eigvec,eigval] = truncevd(symmat)
% truncpartialevd Calculate the first n eigenvalues for a symmetric matrix,
% then truncate to positive eigenvalues
%INPUT:
% symmat    a symmetric matrix
% n         number of eigenvalue/eigenvectors pairs to extract
%OUTPUT: 
% eigvec    eigenvectors associated with positive eigenvalues
% eigval    diagonal matrix of positive eigenvalues
[vec,val] = eig (symmat);
val = diag(val);

[~,sindex] = sort(val,'descend');
vec = vec(:,sindex);
val = val(sindex);
index = val> (0.01*val(1));
eigvec = vec(:,index);
eigval = diag(val(index));
end
