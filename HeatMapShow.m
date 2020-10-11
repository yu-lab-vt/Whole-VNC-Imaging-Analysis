function [upperBound,monMap] = HeatMapShow(zAvg, colorLUTHeatMap,upperBoundInput,cutFlag)

if nargin<4
    cutFlag = 1;
end
[h,w,z] = size(zAvg);
zAvg = zAvg(:,:,z:-1:1);
if z==13
    blks = 4;
else
    if cutFlag
        zAvg(:,:,[1:15,end-14:end]) = [];
        blks = [6 11];%11==> number of rows, 6==> number of columns
    else
        blks = [10 10];
    end
end
monMap = zeros(blks(2)*h,blks(1)*w);
for i = 1:size(zAvg,3)
    zStack = i;
    ovBiIm = zAvg(:,:,zStack);
    [yh,xw] = ind2sub([blks(1),blks(2)],i);
    monMap(1+(xw-1)*h:xw*h, 1+(yh-1)*w:yh*w) = ovBiIm;
end

if nargin>2
    upperBound = upperBoundInput;
else
    tmpMonMap = monMap(:);
    smon = sort(tmpMonMap(~isnan(tmpMonMap)),'descend');
    upperBound = smon(round(0.01*length(smon)));
end
figure;imagesc(monMap);hold on;
if nargin>1
    colormap(colorLUTHeatMap(:,2:4)/255);
end
caxis([0 upperBound]);hold on;
set(gca,'YTick',[]);hold on;
set(gca,'XTick',[]);hold on;
pbaspect([4 3 3]);hold off;