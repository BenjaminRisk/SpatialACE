function [out] = allvertexmwle_ACE(mleresults,h,lat,long,familyst,useparfor)

% mleresults: output allvertexmle_ACE
% familyst: 
    %MZtp1 
    %MZtp2 
    %DZtp1 
    %DZtp2 
    %MDti 
if nargin<6
    useparfor = true;
end
% TO DO: Add MWLE UNDER NULL THAT GENETIC EFFECTS EQUAL ZERO and INFERENCE
R = mleresults.fixefresid;
[N,nVertex] = size(R);
sigmasqA = zeros(1,nVertex);
sigmasqC = zeros(1,nVertex);
sigmasqE = zeros(1,nVertex);
eFlag_a = zeros(1,nVertex);

mlesigmasqA = mleresults.sigmasqA;
mlesigmasqC = mleresults.sigmasqC;
mlesigmasqE = mleresults.sigmasqE;
if useparfor
    parfor v=1:nVertex
    %for v=1:nVertex    
        initvalues = [mlesigmasqA(v),mlesigmasqC(v),mlesigmasqE(v)];
        loginitvalues = log(initvalues);
        %dist = double(distance(lat(v),long(v),lat,long));   
        dist = double(mygreatcirc(lat(v),long(v),lat,long));   
        results = mwle_ACE(R,dist,h,familyst,loginitvalues);
        sigmasqA(v) = results.sigmasqA;
        sigmasqC(v) = results.sigmasqC;
        sigmasqE(v) = results.sigmasqE;

        eFlag_a(v) = results.exitflag_a;
    end

else
    for v=1:nVertex    
    initvalues = [mlesigmasqA(v),mlesigmasqC(v),mlesigmasqE(v)];
    loginitvalues = log(initvalues);
    %dist = double(distance(lat(v),long(v),lat,long));   
    dist = double(mygreatcirc(lat(v),long(v),lat,long));   
    results = mwle_ACE(R,dist,h,familyst,loginitvalues);
    sigmasqA(v) = results.sigmasqA;
    sigmasqC(v) = results.sigmasqC;
    sigmasqE(v) = results.sigmasqE;
    
    eFlag_a(v) = results.exitflag_a;
    end

end

% annoying programming for parfor compatibility
out.sigmasqA = sigmasqA;
out.sigmasqC = sigmasqC;
out.sigmasqE = sigmasqE;

out.eFlag_a = eFlag_a;

out.h2 = (sigmasqA./(sigmasqA+sigmasqC+sigmasqE))';
end
