function Map2dTo3d = FindPositionOfActData(vidIn, ImsIn, mingap)
% We have two 3D data, one is activity data "vidIn" with ~10 stacks, the 
% other is anatomical data "ImsIn" with ~100 stacks
% We are aligning 13 stacks to 100 stacks to find the relative accurate
% depth of each stack.


if nargin<3
    mingap = 1;
end
if size(ImsIn,4)>1 % anatomical data contains multiple frames
    % the last time point is blank
	Ims = double(ImsIn(:,:,1:size(ImsIn,3)-1,:));% use the first t-1 time points to do registration
else
    [h,w,z] = size(ImsIn);
    Ims = reshape(ImsIn,[h,w,1,z]);
end
% do a smooth to the data
Ims = squeeze(mean(Ims,3));

Ims = smooth3(Ims, 'gaussian', [3,3,3]);
[h,w,st3d] = size(Ims);
tp3d = 1;
Ims = reshape(Ims, h,w,tp3d, st3d);


%liffilepath = [rtpath, vid2D];
%vidIn = imreadBF(liffilepath, [], 1:10,[], 2);% we only need to first 10 time points for alignment
vid = squeeze(mean(double(vidIn(:,:,1:min(10,size(vidIn,3)),:)),3));
[~,~,st2d] = size(vid);
%vid = scale_image(vid,0,1);
%Ims = scale_image(Ims,0,1);

maxProjIms = max(Ims,[],4);
maxProjVid = squeeze(max(vid,[],3));

par.transform = 'euclidean';
par.levels = 2;
par.iterations = 50; %iterations per level
ECCWarp = cell(tp3d,1);
for i=1:tp3d
    IRef = maxProjIms(:,:,i);
    ECCWarp{i} = iat_ecc(maxProjVid, IRef, par);
end

Map2dTo3d = nan(st2d,1);
for j = 1:st2d
    distVec = nan(st3d,1);
    if j>1
        startStack = min(Map2dTo3d(j-1)+mingap,st3d);
    else
        startStack = 1;
    end
    distVec(startStack:end) = 0;
    for i = 1:tp3d
        [wimECC, ~] = iat_inverse_warping(vid(:,:,j), ECCWarp{i}, 'euclidean', 1:w, 1:h);
        for k = startStack:st3d
            distVec(k) = distVec(k)+norm(wimECC-Ims(:,:,i,k),'fro');
        end
    end
    [~, Map2dTo3d(j)] = nanmin(distVec);
    %figure;plot(distVec)
end
