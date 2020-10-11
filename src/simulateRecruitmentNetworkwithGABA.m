function [res, ttest_pvalue] = simulateRecruitmentNetworkwithGABA(...
    numCell, ... % number of neurons in the simulation
    inSti, ...   % input stimulation levels
    saturateIn, ... % saturated stimulation level that a neuron responses
    minimalIn, ...  % minimum stimulation level that a neuron responses
    saturateResponse, ... % saturated response level
    nslevels, ...   % noise level
    numExp ...      % number of experiments
    )
% A function to simulate the recruitment network.
if nargin == 0
    numCell = 58;
    
    inSti = 0:0.1:100;%20:0.1:80;%[0.1:0.2:15];%, 2:20];%[0:0.1:1, 1:10:100];, 11:2:20% stimulation level
    saturateIn = [50,80];
    minimalIn = [20,50];
    
    saturateResponse = 0.5; % maximal calcium dff level ==> 50% dff
    
    numExp = 1000;
    nslevels = 0.2; % can set multiple noise level for testing
end
roc_case = [];
roc_ctr = [];
data_cv = [];
gt_cv = [];
valid_range = cell(length(nslevels), 2);
for nn = 1:length(nslevels)
    nslevel = nslevels(nn);
    
    outCell = cell(numExp,2);
    for j=1:numExp
        curSatuarteIn = saturateIn(1)+(saturateIn(2)-saturateIn(1))*rand(numCell,1);% 5-15 saturate intensities
        curMinimalIn = minimalIn(1)+(minimalIn(2)-minimalIn(1))*rand(numCell,1);
        % GABA inhibition
        GABAThres = curMinimalIn + (curSatuarteIn-curMinimalIn).*rand(numCell,1)*(minimalIn(1)/inSti(end));
        GABAInhibit = saturateResponse.*((GABAThres-curMinimalIn)./(curSatuarteIn-curMinimalIn));
        GABAInhibit = sum(GABAInhibit);
        outResponse = zeros(8, length(inSti));
        % number increas with response increase
        positive_thres = saturateResponse*numCell/2;
        noiseStandard = saturateResponse*(median(inSti)-mean(minimalIn))./(mean(saturateIn)-mean(minimalIn));
        %noiseStandard = saturateResponse*(median(curSti)-curMinimalIn)./(curSatuarteIn-curMinimalIn);
        for i=1:length(inSti)
            curSti = inSti(i);
            curResponse = saturateResponse*...
                (max(curSti-curMinimalIn, 0)./(curSatuarteIn-curMinimalIn));
            curResponse(curResponse>saturateResponse) = saturateResponse;
            
            gt_res = sum(curResponse);
            curResponse = curResponse + noiseStandard.*randn(numCell,1)*nslevel; % noise #1
            
            curResponse = sum(curResponse);
            if curResponse<GABAInhibit
                curResponse = 0;
            else
                curResponse = curResponse - GABAInhibit;
            end
            if gt_res<GABAInhibit
                gt_res = 0;
            else
                gt_res = gt_res - GABAInhibit;
            end
            curResponse = curResponse + positive_thres*randn(1,1)*nslevel; % noise #2
            outResponse(1, i) = curResponse;
            outResponse(2, i) = gt_res;
            if curResponse>=positive_thres
                outResponse(3, i) = 1;
            else
                outResponse(3, i) = 0;
            end
            if gt_res>=positive_thres
                outResponse(4, i) = 1;
            else
                outResponse(4, i) = 0;
            end
            
        end
        
        % no cell number increase, no GABA
        positive_thres = saturateResponse*numCell/2;
        for i=1:length(inSti)
            curSti = inSti(i);
            curResponse = ones(numCell,1).*saturateResponse.*...
                (max(curSti-minimalIn(1),0)./(saturateIn(2)-minimalIn(1)));
            curResponse(curResponse>saturateResponse) = saturateResponse;
            
            gt_res = sum(curResponse);
            curResponse = curResponse + noiseStandard.*randn(numCell,1)*nslevel; % noise #1
            curResponse = sum(curResponse);
            curResponse = curResponse + positive_thres*randn(1,1)*nslevel; % noise #2
            
            outResponse(5, i) = sum(curResponse);
            outResponse(6, i) = gt_res;
            if curResponse>=positive_thres
                outResponse(7, i) = 1;
            else
                outResponse(7, i) = 0;
            end
            if gt_res>=positive_thres
                outResponse(8, i) = 1;
            else
                outResponse(8, i) = 0;
            end
            
        end
        % display
        %         f = figure('Visible','off');
        %         plot(outResponse(1, :)); hold on;
        %         plot(outResponse(2, :)); hold on;
        %         plot(outResponse(5, :)); hold on;
        %         plot(outResponse(6, :)); hold on;
        %         %plot(1:length(outResponse(6, :)), zeros(length(outResponse(6, :)),1)+positive_thres);hold off;
        %         legend('case','case_gt', 'ctr', 'ctr_gt');
        %         title(['response curve under noise std ', num2str(nslevel)]);
        %         saveas(f, ['response curve under noise std ', num2str(nslevel), '.tif'], 'tif');
        outCell{j,1} = outResponse(1:4,:);
        outCell{j,2} = outResponse(5:8,:);
        %allResponseLevel{3} = cat(2,allResponseLevel{3}, outResponse(:,3));
        %% generate the curves, by setting numExp = 1
        data_cv = cat(1, data_cv, outResponse(1, :), outResponse(5, :));
        gt_cv = cat(1, gt_cv, outResponse(2, :), outResponse(6, :));
    end
    
    allResponseLevel{1} = cat(1,outCell{:,1});
    allResponseLevel{2} = cat(1,outCell{:,2});

    
    %% find the valid range
    valid_sti = 1:numel(inSti);
    %% way 1: use start point and end point
    st_pt = 350;
    ed_pt = 650;
    
    %% start calcultae fpr and tpr
    disp([st_pt, ed_pt]);
    acc_range = 0;%[0.1:0.1:0.9];
    for i=1:2
        tmp = allResponseLevel{i}(3:4:end,valid_sti);
        tmpgt = allResponseLevel{i}(4:4:end,valid_sti);
        valid_range{nn, i} = find_valid_range(tmp, tmpgt, acc_range, inSti);
        
        tmp = tmp(:,st_pt:ed_pt);
        tmpgt = tmpgt(:,st_pt:ed_pt);
        tpr = length(find(tmp==1 & tmpgt==1))/length(find(tmpgt == 1));
        fpr = length(find(tmp==1 & tmpgt==0))/length(find(tmpgt == 0));
        
        acc = length(find(tmp==tmpgt))/numel(tmpgt);
        %disp([tpr, fpr]);
        if i==1
            roc_case = cat(1, roc_case, [fpr tpr, acc]);
        else
            roc_ctr = cat(1, roc_ctr, [fpr tpr, acc]);
        end
    end
end
res = cat(2, roc_case, roc_ctr);
% figure;plot(roc_case(:,1), roc_case(:,2)); hold on;
% plot(roc_ctr(:,1), roc_ctr(:,2)); legend('case', 'control');

avg_rg_ctr = zeros(size(valid_range,1), 2*numel(valid_range{1}));
avg_rg_case = zeros(size(valid_range,1), 2*numel(valid_range{1}));

ttest_pvalue = zeros(size(valid_range,1), numel(valid_range{1}));

for i=1:size(valid_range,1)
    for j=1:numel(valid_range{1})
        avg_rg_case(i,2*j-1) = mean(valid_range{i,1}{j});
        avg_rg_case(i,2*j) = std(valid_range{i,1}{j})...
            /sqrt(length(valid_range{i,1}{j}));
        
        avg_rg_ctr(i,2*j-1) = mean(valid_range{i,2}{j});
        avg_rg_ctr(i,2*j) = std(valid_range{i,2}{j})...
            /sqrt(length(valid_range{i,2}{j}));
        [~,p] = ttest(valid_range{i,1}{j}, valid_range{i,2}{j});
        ttest_pvalue(i,j) = p;
    end
end