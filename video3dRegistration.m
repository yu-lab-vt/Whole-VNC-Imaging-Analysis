function out3dVideo = video3dRegistration(in3dVideo)
% registrate a 3dvideo, in x, y direction
% input in3dVideo: a 4D matrix, h*w*t*z
% way one: registrate the max projection and apply to all stacks
% way two: registrate all stacks independently
in3dVideo = double(in3dVideo);
in2dVideo = max(in3dVideo, [], 4);
in2dVideo = squeeze(in2dVideo);
ECCWarp = video2dRegistration(in2dVideo);
[h,w,t,z] = size(in3dVideo);
out3dVideo = zeros(h,w,t,z,'like',in3dVideo);
for j = 1:z
    out3dVideo(:,:,1,j) = in3dVideo(:,:,1,j);
end
parfor i = 2:t
    for j=1:z
        img = in3dVideo(:,:,i,j);
        [wimECC, ~] = iat_inverse_warping(img, ECCWarp{i}, 'euclidean', 1:w, 1:h);
        out3dVideo(:,:,i,j) = wimECC;
    end
end