function alignedResMap = VNC_response_map_generating(...
    anatomical_data,... % 3d anatomical data (h*w*z)
    activity_data, ...  % 4d activity data (h*w*t*z)
    stimulation_t, ...  % the time of stimulation (frame number)
    template_anatomical_data, ... % 3d template anatomical data
    phase_one_length, .... % length of time for measure response
    fg_threshold ...    % a loose threshold for foreground detection
    )
% This function measures the response level of a larva's VNC to stimulus.
% The response level is represented by z-score. The output data is aligned
% to the template and thus comparable across different samples.

% All the anatomical data and activity data is assumed to be uint16. If not,
% you may want to change fg_threshold.

% contact: ccwang AT vt DOT edu

if nargin == 5
    fg_threshold = 500;
end
addpath(genpath('registration_tools'));
anatomical_data = single(anatomical_data);
activity_data = single(activity_data);
template_anatomical_data = single(template_anatomical_data);
%% Align the anatomical data to template
fprintf('Align the anatomical data to template...\n');
% second, align anatomical data to template with b-spline
Options = [];
Options.MaskMoving = anatomical_data>fg_threshold;
Options.MaskStatic = template_anatomical_data>fg_threshold;
% The masks in Options is not essential. For most case, the registration
% results show no clear difference without masks. The reason we have it
% here is because that the b-spline algorithm is not stable without masks.
[Ireg,Grid,Spacing,M,B,F] = image_registration(anatomical_data,...
    template_anatomical_data, Options);
% save warping parameters
regOut = cell(6,1);
regOut{1} = Ireg;
regOut{2} = Grid;
regOut{3} = Spacing;
regOut{4} = M;
regOut{5} = B;
regOut{6} = F;

%% Align the acitivty data to its own anatomical data
fprintf('Align the activity data to anatomical data...\n');
% first, correct the drifting of the acitvity data
activity_data = video3dRegistration(activity_data);
% second, aligth to its anatomical data
Map2dTo3d = FindPositionOfActData(activity_data, anatomical_data, ...
    phase_one_length);

%% Generate the response map
% using interpolation and the warping function from second and third steps,
% we can upsample the response map to the same size of template
% registrate the sparse 3d data by registrating maxprojected 2d data
fprintf('Get the response level (z-score) of voxels from the activity data...\n');
[h,w,~,z] = size(activity_data);
bgMap = zeros(h,w);
validMap = cell(z,1);
for j=1:z
    ImStack = squeeze(activity_data(:,:,:,j));
    validMap{j} = ...
        mean(ImStack(:,:,max(stimulation_t-10,1):stimulation_t-1),3)...
        >=fg_threshold;
    bgMap(validMap{j}) = 1;% foreground
end
resMap = zeros(h,w,z);
corrMap = zeros(h,w,z);
difIm = cell(z,1);
shiftMap = cell(z,1);% rise ratio in each seconds
for j=1:z
    ImStack = squeeze(activity_data(:,:,:,j));
    [resMap(:,:,j), corrMap(:,:,j), difIm{j}, ...
        shiftMap{j}] = detectImmediateResV2(ImStack, stimulation_t, ...
        validMap{j},3);
end
ResponseInfo{1} = resMap;
ResponseInfo{2} = corrMap;
ResponseInfo{3} = difIm;
ResponseInfo{4} = shiftMap;

%% Align the response to template
% use the backward transformation with interpolation
zMap = ResponseInfo{1};
% target coordiate
movMap = regOut{5};
[tmpH,tmpW,tmpZ] = size(template_anatomical_data);
[txx,tyy,tzz] = meshgrid(1:tmpW,1:tmpH,1:tmpZ);
txx = txx+squeeze(movMap(:,:,:,2));% !!! attention on the order
tyy = tyy+squeeze(movMap(:,:,:,1));
tzz = tzz+squeeze(movMap(:,:,:,3));
% source coordiate
[srcH,srcW,~] = size(zMap);
[sxx,syy] = meshgrid(1:srcW,1:srcH);
XX = repmat(sxx,1,1,length(Map2dTo3d));
YY = repmat(syy,1,1,length(Map2dTo3d));
ZZ = reshape(Map2dTo3d(:),1,1,[]);
ZZ = repmat(ZZ,srcH,srcW,1,1);
alignedResMap = interp3(XX,YY,ZZ,zMap, txx,tyy,tzz, 'linear');

