function corrMap = getArResCorr2( tsts, validMap, cBias )
%GETRESCORR Compute active region residual correlation
% support NaN values
% credit: Yizhi Wang
if ~exist('cBias','var')
    cBias = 0;
end

if size(tsts,3)>1
    [~,~,nTps] = size(tsts);
    tsts = reshape(tsts,[],nTps);
    valIdx = find(validMap);
    tsts = tsts(valIdx,:);
end
[nPix,nTps] = size(tsts);
validIdxMap = zeros(size(validMap));
validIdxMap(validMap>0) = 1:sum(validMap(:)>0);
[H,W] = size(validMap);

% change 1
% 8-neighbor
%     hofst = [-1,0,1,-1,1,-1,0,1];
%     wofst = [-1,-1,-1,0,0,1,1,1];
% 16 neighbor 5*5-3*3
wofst = [-2,-2,-2,-2,-2, 2, 2, 2, 2, 2,-1,-1, 0, 0, 1, 1];
hofst = [-2,-1, 0 1, 2,-2,-1, 0, 1, 2,-2, 2,-2, 2,-2, 2];
numNeib = length(hofst);
[iy,ix] = find(validMap);

ny = repmat(iy,1,numNeib);
ny = bsxfun(@plus,ny,hofst);
nx = repmat(ix,1,numNeib);
nx = bsxfun(@plus,nx,wofst);

Nmapidx = (nx-1)*H+ny;
Invalididx = ny<1 | ny>H | nx<1 | nx>W;
NullIdx = find(validMap==0,1);
if isempty(NullIdx)
    %NullIdx = H*W+1;
    Nmapidx(Invalididx) = 0;% change outliers to null idx
    Nmapidx = Nmapidx';
    Nmapidx = Nmapidx(:);
    Nidx = zeros(length(Nmapidx),1);
    valNmapidx = Nmapidx>0;
    Nidx(valNmapidx) = validIdxMap(Nmapidx(valNmapidx)); % index in dif0
    Nidx(~valNmapidx) = 0;
else
    Nmapidx(Invalididx) = NullIdx;% change outliers to null idx
    Nmapidx = Nmapidx';
    Nidx = validIdxMap(Nmapidx(:)); % index in dif0
end
Nidx(Nidx==0) = nPix+1;% set all invalid idx to nPix+1
tmpdif0 = cat(1, tsts,zeros(1,nTps));% the last one is all zeros
Ncurve = tmpdif0(Nidx,:);
Ncurve = reshape(Ncurve,numNeib,length(iy),nTps);

% way 1
if 0
    Ndif0 = squeeze(sum(Ncurve,1));
    if size(Ndif0,1)==1 || size(Ndif0,2)==1
        Ndif0 = Ndif0(:)';
    end
    Ndif0demean = bsxfun(@minus,Ndif0,mean(Ndif0,2));
    dif0demean = bsxfun(@minus,tsts,mean(tsts,2));
    c1 = sum(Ndif0demean.*dif0demean,2)./...
        (sqrt(var(Ndif0,[],2).*var(tsts,[],2)));
    c1 = c1/(nTps-1);
    [scc1, scpos] = sort(c1);
    if scc1(1) < -0.99 % this one maybe reference pixel
        c1(scpos(1)) = scc1(2);
    end
else% way 2
    dif0demean = bsxfun(@minus,tsts,mean(tsts,2));
    varTst = var(tsts,[],2);
    c1 = zeros(size(dif0demean,1),1);
    parfor i=1:size(Ncurve,1)
        Ndif0Single = Ncurve(i,:,:);
        Ndif0Single = squeeze(Ndif0Single);
        Ndif0demean = bsxfun(@minus,Ndif0Single,mean(Ndif0Single,2));
        c1Single = sum(Ndif0demean.*dif0demean,2)./...
            (sqrt(var(Ndif0Single,[],2).*varTst));
        c1Single = c1Single/(nTps-1);
        [scc1, scpos] = min(c1Single);
        if scc1 < -0.99 % this one maybe reference pixel
            c1Single(scpos) = 1;
            c1Single(scpos) = min(c1Single);
        end
        c1 = c1+c1Single;
    end
    c1 = c1/size(Ncurve,1);
end
    
if cBias~=0
    c1 = max(c1-cBias,0);  % bias for background correlation
end
corrMap = validMap*0;
corrMap(iy+(ix-1)*H) = c1;
    
    % % z1 = 0.5*log((1+c1)./(1-c1))*sqrt(nTps-3);
    % % if nargout<3
    % %     p1 = 0;
    % % else
    % %     p1 = 1-normcdf(z1,0,1);
    % % end
end





