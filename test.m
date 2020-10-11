%% test function for a single VNC data analysis
% download data for testing
if exist('test_data/template_anatomical_data.mat','file') && ...
    exist('test_data/test_sample_data.mat','file')
else
    fprintf('Please download the testing data and template first.\n');
    fprintf('The files can be found here:\n');
    url = 'https://drive.google.com/drive/folders/1w13FSuk6rYh07wa9RLszVn4KzqJUFIng?usp=sharing';
    fprintf('%s\n', url);
    fprintf('After downloading, put them in the folder of ''test_data''.\n');
    return;
end
% load template data
load('test_data/template_anatomical_data.mat');
% load sample data
load('test_data/test_sample_data.mat');

% start test
phase_one_length = 5; 
alignedResMap = VNC_response_map_generating(test_sample_anatomical_data, ...
    test_sample_activity_data, stimulus_frame, template_anatomical, ...
    phase_one_length);
% display results
HeatMapShow(alignedResMap, LUTHeatMap);


