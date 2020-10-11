function [resMap, corrMap, difIm, shiftMap] = ...
    detectImmediateResV2(ImStack,stiTime,validMap, numTps)
% map response level to z-score. We use order statistics to correct the
% possible bias caused by the max operation in the phase one.

if nargin == 3
    numTps = 5;
end

[h,w, t] = size(ImStack);
[ImStackVST, varEst ] = estNoiseStbVideoSimple(ImStack,validMap);
ImStackVST = smoothMovie( ImStackVST,[],2,2,0);% smooth method 2, sigma=1 ImStack;%
varEst = varEst*0.0218; % we use gaussian filter(fspecial('gaussian',9,2);) to smooth, its square sum is 0.0218

sigEst = sqrt(varEst);
corrMap = getArResCorr2( ImStack(:,:,max(stiTime-1,1):min(stiTime+20, t)), ones(h,w));

charCurves = reshape(ImStackVST,[],t);
dff = zeros(size(charCurves));
f0 = zeros(h*w,1);
% a sliding window based method to determine f0
parfor i=1:h*w
    curCv = charCurves(i,:);
    [~,minvarPos] = min(movmean(curCv-min(curCv),5).*movvar(curCv,5));
    stPos = max(minvarPos-2,1);
    endPos = min(minvarPos+2,t);
    f0(i) = mean(charCurves(i,stPos:endPos));
    dff(i,:) = (charCurves(i,:)-f0(i))/f0(i);
end

difIm = single(reshape(dff, h,w,t));
difIm = difIm(:,:,stiTime:min(stiTime+10,t));% previous is stiTime:stiTime+10
shiftMap = single((ImStackVST(:,:,stiTime:min(stiTime+10,t))-ImStackVST(:,:, stiTime-1:min(stiTime+10,t)-1))./(sqrt(2)*sigEst));

resVec = zeros(h*w,1);
difVec = zeros(h*w,1);

if stiTime+numTps>t
    numTps = t-stiTime;
end
if numTps ~= 5 % use simulation to get the parameter of order statistics
    odSt = randn(100000,numTps);
    odSt = max(odSt,[],2);
    vM = var(odSt);
    mM = mean(odSt);
else
    vM = 0.447534;
    mM = 1.162964;
end

for i=1:h*w
    estimated_sig = sigEst/f0(i);
    
    difVec(i) = max(dff(i,stiTime+1:stiTime+numTps))-dff(i,stiTime);
    mu = mM*estimated_sig;
    sigma = sqrt(1+vM)*estimated_sig;
    resVec(i) = (difVec(i)-mu)/sigma;
    
end

resMap = single(reshape(resVec,h,w));
end