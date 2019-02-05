function [out] = smw(K,T,S)
%Sherman-morrison woodbury formula for
%the inverse of (KTK'+S) adapted to PSD T
invS = inv(S);
[evec,evalue] = truncevd(T);
K = K*evec;
out = invS - invS*K*inv(inv(evalue)+K'*invS*K)*K'*invS;
end

