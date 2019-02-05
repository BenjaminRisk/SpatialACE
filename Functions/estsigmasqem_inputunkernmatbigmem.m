function [smsigmasqem] = estsigmasqem_inputunkernmatbigmem(R,sigmasqA,sigmasqC,sigmasqE,unkernmat)
%ESTSIGMAEM 
%       Estimates the covariance matrix for pooled additive genetic, 
%       common environmental, unique environmental.
%
%INPUT:
% R: N x V matrix; residual matrix
% sigmasqA
% sigmasqC
% sigmasqE: estimated MEASUREMENT ERROR PLUS UNIQUE ENVIRONMENTAL VARIANCES
% unkernmat

%OUTPUT:
% sigmasqem

% create raw covariance estimates:
[N,nVertex] = size(R);
if N>nVertex 
    warning('N > V -- check that R is N x V')
end


S0 = R'*R./N;
S0 = S0 - diag(diag(S0));
smSaceg = zeros(nVertex,1); 
rowsums = sum(unkernmat);
bigW = rowsums*rowsums' - unkernmat*unkernmat;

smSaceg = diag(unkernmat*S0*unkernmat)./diag(bigW);

nna = sum(isnan(smSaceg));
if nna
    warning(['bandwidth contains no neighbors for ' num2str(nna) ' vertices'])
end

inputSace = sigmasqA+sigmasqC+sigmasqE;
smsigmasqem = inputSace - smSaceg';
end


