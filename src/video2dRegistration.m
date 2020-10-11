function [ECCWarp, Iout] = video2dRegistration(in2dVideo)
% registrate a video
[h,w,t] = size(in2dVideo);
%% using IAT to do registration
IRef = in2dVideo(:,:,1);
ECCWarp = cell(t,1);
par.transform = 'euclidean';
par.levels = 2;
par.iterations = 50; %iterations per level

parfor i=2:t
    img = in2dVideo(:,:,i);
    ECCWarp{i} = iat_ecc(img, IRef, par);
end
if nargout>1
    Iout = zeros(h,w,t);
    Iout(:,:,1) = IRef;
    [M,N] = size(IRef);
    parfor i=2:t
        img = in2dVideo(:,:,i);
        [wimECC, ~] = iat_inverse_warping(img, ECCWarp{i}, par.transform, 1:N, 1:M);
        %[~, grayerrorECC] = iat_error2gray(IRef,wimECC,suportECC);
        %figure;imshowpair(uint8(wimECC), IRef,'Scaling','joint');
        Iout(:,:,i) = wimECC;
    end
end