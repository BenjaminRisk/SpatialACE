function [newsigmasqem] = impute_missing_sigmasqem(sigmasqem,lat,long)
%IMPUTE_MISSING_SIGMASQEM
%       Imputes measurement error for missing values with the values of the
%       closest non-missing vertex
%INPUT:
% sigmasqem
% lat
% long

%OUTPUT:
% sigmasqem with no missings

newsigmasqem = sigmasqem;
indices = find(isnan(sigmasqem));
    for v=1:length(indices)
        dists = mygreatcirc(lat(indices(v)),long(indices(v)),lat,long);
        [~,sindices] = sort(dists);
        temp = sigmasqem(sindices);
        a = temp(~isnan(temp));
        newsigmasqem(indices(v)) = a(1);
    end
end


