function [ vidVST, var1, pEst0 ] = estNoiseStbVideoSimple( vid,nMap,tps )
%ESTNOISESTBVIDEO Estimate the noise parameters and stabilize the movie
% credit: Yizhi Wang
if ~exist('tps','var')
    tps = 1:size(vid,3);
end

% noise estimation
nTps = size(vid,3);
vidVec = reshape(vid,[],nTps);
idx = nMap==0;
vidVecSel = vidVec(idx,tps);

xMean = mean(vidVecSel,2);
xVar = var(vidVecSel,0,2);

% linear regression
k = xMean\xVar;
b = mean(xVar) - k*mean(xMean);

vidVST = vid*0;

% Stabilization
pEst0 = [k,b];
alpha = pEst0(1);
% sigma2 = max(pEst0(2),2e-3);
sigma2 = max(pEst0(2),1e-6);
t0 = 2/alpha*sqrt(alpha*0+3/8*alpha^2+sigma2);
t1 = 2/alpha*sqrt(alpha*1+3/8*alpha^2+sigma2);
var1 = 1/(t1-t0)^2;
for ii=1:nTps
    dat1 = 2/alpha*sqrt(alpha*vid(:,:,ii)+3/8*alpha^2+sigma2);
    dat1 = (dat1-t0)/(t1-t0);
    vidVST(:,:,ii) = dat1;
end


end

