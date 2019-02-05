function [smsigmasqem] = estsigmasqem_inputunkernmat(R,sigmasqA,sigmasqC,sigmasqE,unkernmat)
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
fswitch = 0;
for v=1:nVertex
  vunkernmat = unkernmat(:,v);
  bigW = rowsums(v)^2 - vunkernmat'*vunkernmat;
  if bigW == 0 
      fswitch = fswitch+1;
      smSaceg(v) = NaN;
  else
    smSaceg(v) = vunkernmat'*S0*vunkernmat./bigW;
  end
end

if fswitch
    warning(['bandwidth contains no neighbors for ' num2str(fswitch) ' vertices'])
end

inputSace = sigmasqA+sigmasqC+sigmasqE;
smsigmasqem = inputSace - smSaceg';
end


