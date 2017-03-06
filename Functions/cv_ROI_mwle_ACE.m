function [out] = cv_ROI_mwle_ACE(mleresults, ROI, hvec, lat, long, familyst,useparfor)
%Benjamin Risk
%6 March 2017
%cv_ROI_mwle_ACE Use cross-validation to select a single bandwidth for all
%vertices specified by the ROI
%Input:
% mleresults: struct output from allvertexmle_ACE
% ROI: V x 1 logical vector indicating which vertices to include in ROI
% hvec: vector of bandwidths to evaluate
% lat: V x 1 vector of latitudes
% long: V x 1 vector of longitudes
% familyst:
    % MZtp1 -- N x 1: ==1 if subject is the first twin for MZ pair
    % MZtp2 -- N x 1: ==1 if subject is the second twin for MZ pair
    % DZtp1 -- N x 1: ==1 if subject is the first twin for DZ pair
    % DZtp2 -- N x 1: ==1 if subject is the second twin for DZ pair
    % MDti -- N x 1: ==1 if subject is a singleton.
    % familyID -- N x 1 : integer-valued vector with the family id for each
    %   subject
% initvalues -- 3 x V : optional initial values. Typically,
% these are the estimates from the unweighted MLE (note they are on  the standard scale, i.e., not logs)
%Output:
% out.hvec: values at which bandwidth evaluated
% out.mse: vector of errors for each hvec

if nargin<7
    useparfor = true;
end

initvalues = [mleresults.sigmasqA;mleresults.sigmasqC;mleresults.sigmasqE];
R = mleresults.fixefresid;

MZtp1 = familyst.MZtp1;
MZtp2 = familyst.MZtp2;
DZtp1 = familyst.DZtp1;
DZtp2 = familyst.DZtp2;
MDti = familyst.MDti;
familyID = familyst.familyID;

[nSubject,nVertex] = size(R);
if nSubject>nVertex 
    warning('N > V -- check that R is N x V')
end
nbw = length(hvec);
nVertexROI = sum(ROI);
K = 5;
mse = zeros(nbw,K,nVertexROI);
famuniq = unique(familyID);
nfam = length(famuniq);
twinsmat = [MZtp1,MZtp2,DZtp1,DZtp2,MDti];

% use cross validation with FAMILY as the unit of observation
cvpart = cvpartition(nfam,'Kfold',K);
vIndex = find(ROI);
if useparfor
    parfor v=1:nVertexROI
    %for v=1:nVertexROI   
        vTemp = vIndex(v);
        tempinit = log(initvalues(:,vTemp));
        %dist = double(distance(lat(vTemp),long(vTemp),lat,long)); 
        dist = double(mygreatcirc(lat(vTemp),long(vTemp),lat,long)); 
        
        for j=1:nbw
            h = hvec(j);
            index = (dist<=h);
            subR = R(:,index);
            wts = biweight(dist(index)/h)/h;
            sumwts = sum(wts);
            for k = 1:K
              % get indices for subjects
              fampart = famuniq(cvpart.training(k));
              subjIndex = ismember(familyID,fampart);
              tTwinsmat = twinsmat(subjIndex,:);
              tSubR = subR(subjIndex,:);
              n1 = sum(tTwinsmat(:,1)); n2 = sum(tTwinsmat(:,3)); n3 = sum(tTwinsmat(:,5));

              SM = tSubR(tTwinsmat(:,1),:).^2 + tSubR(tTwinsmat(:,2),:).^2; XM = 2*tSubR(tTwinsmat(:,1),:).*tSubR(tTwinsmat(:,2),:); % n1 x M
              SD = tSubR(tTwinsmat(:,3),:).^2 + tSubR(tTwinsmat(:,4),:).^2; XD = 2*tSubR(tTwinsmat(:,3),:).*tSubR(tTwinsmat(:,4),:); % n2 x M
              SMD = tSubR(tTwinsmat(:,5),:).^2; % n3 x M

              sumSM = sum(SM)*wts;
              sumXM = sum(XM)*wts;
              sumSD = sum(SD)*wts;
              sumXD = sum(XD)*wts;
              sumSMD = sum(SMD)*wts;

              options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton');
              [paramsA] = fminunc(@(params) wtdloglik_ACE(params, sumSM, sumXM, sumSD, sumXD, sumSMD, sumwts, n1, n2, n3), ...
                   tempinit, options);     
              trainIndex =  logical(1 - subjIndex);
              tempR = R(trainIndex,vTemp);
              a = sum(exp(paramsA));
              mse(j,k,v) = mean((tempR.^2 - a).^2);
            end
        end
    end
else
    for v=1:nVertexROI   
        vTemp = vIndex(v);
        tempinit = log(initvalues(:,vTemp));
        %dist = double(distance(lat(vTemp),long(vTemp),lat,long));   
        dist = double(mygreatcirc(lat(vTemp),long(vTemp),lat,long));   
        
        for j=1:nbw
            h = hvec(j);
            index = (dist<=h);
            subR = R(:,index);
            wts = biweight(dist(index)/h)/h;
            sumwts = sum(wts);
            for k = 1:K
              % get indices for subjects
              fampart = famuniq(cvpart.training(k));
              subjIndex = ismember(familyID,fampart);
              tTwinsmat = twinsmat(subjIndex,:);
              tSubR = subR(subjIndex,:);
              n1 = sum(tTwinsmat(:,1)); n2 = sum(tTwinsmat(:,3)); n3 = sum(tTwinsmat(:,5));

              SM = tSubR(tTwinsmat(:,1),:).^2 + tSubR(tTwinsmat(:,2),:).^2; XM = 2*tSubR(tTwinsmat(:,1),:).*tSubR(tTwinsmat(:,2),:); % n1 x M
              SD = tSubR(tTwinsmat(:,3),:).^2 + tSubR(tTwinsmat(:,4),:).^2; XD = 2*tSubR(tTwinsmat(:,3),:).*tSubR(tTwinsmat(:,4),:); % n2 x M
              SMD = tSubR(tTwinsmat(:,5),:).^2; % n3 x M

              sumSM = sum(SM)*wts;
              sumXM = sum(XM)*wts;
              sumSD = sum(SD)*wts;
              sumXD = sum(XD)*wts;
              sumSMD = sum(SMD)*wts;

              options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton');
              [paramsA] = fminunc(@(params) wtdloglik_ACE(params, sumSM, sumXM, sumSD, sumXD, sumSMD, sumwts, n1, n2, n3), ...
                   tempinit, options);     
              trainIndex =  logical(1 - subjIndex);
              tempR = R(trainIndex,vTemp);
              a = sum(exp(paramsA));
              mse(j,k,v) = mean((tempR.^2 - a).^2);
            end
        end
    end
end

out.hvec = hvec;
mseV = sum(mse,3);
mseV = sum(mseV,2);

out.mse = mseV;

[~,b] = min(out.mse);
out.hvecmin = hvec(b);

end

