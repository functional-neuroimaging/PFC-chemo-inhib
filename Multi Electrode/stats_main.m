% routine to check the results for spectrum and coherency, requires the
% results from the main routine, related to figure 5 and supplementary
% figure 9

clear
clc
close all

%% 
Mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a','fr200209a'...
    ,'fr200327a','fr200330a','fr200331a','fr200331b'}; % name of the recordings
traces = {21:30,19:28,31:40,31:40,31:40,31:40,31:40,31:40,31:40}; % traces you want to analyze
results_folder = '';% put the directory where the results are saved
ROI = {'PFC', 'Rs', 'MDTh'}; % regions recorded
n_exp = 5; % the first 5 are exp
%% results
S_pool_one = cell(1,length(Mouse_names));
C_pool_one = S_pool_one;
fc_all = C_pool_one;
fs_all = S_pool_one;

for MN=1:length(Mouse_names)
    mouse_name = Mouse_names{MN};
    path = strcat(results_folder,'\',mouse_name);
    n_traces = traces{MN};

    for i=1:length(ROI)
        S_all.(ROI{i}) = [];
        S_pool_one{MN}.(ROI{i}) = [];
    end
    C_pool_one{MN} = cell(3,1);
    t_all_s = 0;
    C_all = cell(1,3);
    
    for trace=n_traces(1):n_traces(end)
        file_name = strcat(path,'\',mouse_name,'_',num2str(trace),'_results','.mat');
        Data = load(file_name);
        for i=1:length(ROI)
            S_all.(ROI{i}) = [S_all.(ROI{i}),Data.S_normalized.(ROI{i}).CNO];
            alpha_s = Data.S_normalized.(ROI{i}).CNO;
            S_pool_one{MN}.(ROI{i}) = [S_pool_one{MN}.(ROI{i}),mean(alpha_s,2)];
            C_all{i} = [C_all{i};Data.C_normalized{i}];
            alpha_c = Data.C_normalized{i};
            C_pool_one{MN}{i} = [C_pool_one{MN}{i};alpha_c];
        end
        t_all_s = [t_all_s;Data.t_s+t_all_s(end)];
    end
    t_all_s(1) = [];
    t_all_s = t_all_s./60;
    tmp_lag_S = [];
    tmp_lag_C = [];
    for i=1:length(ROI)
        tmp = S_all.(ROI{i});
        tmp_c = C_all{i}';
    end
    fc_all{MN} = Data.f;
    fs_all{MN} = Data.f_s;
end
titles_c = Data.titles_C;
bands = {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70]};
labels = {'infraslow','Delta','Theta','Alpha','Beta','LG'};
%% statistics
for j=1:length(ROI)
    for band=1:length(bands)
        S_exp_to_test = [];
        S_ctrl_to_test = [];
        C_exp_to_test = [];
        C_ctrl_to_test = [];
        for i=1:length(Mouse_names)
            tmp_S = S_pool_one{i}.(ROI{j})(fs_all{MN} > bands{band}(1) & fs_all{MN} <= bands{band}(2),:);
            tmp_S = median(tmp_S);
            tmp_C = C_pool_one{i}{j}(:,fc_all{MN} > bands{band}(1) & fc_all{MN} <= bands{band}(2));
            tmp_C = median(tmp_C,2);
            tmp_C = tmp_C';
            if i<=n_exp
                S_exp_to_test = [S_exp_to_test,tmp_S];
                C_exp_to_test = [C_exp_to_test,tmp_C];
            else
                S_ctrl_to_test = [S_ctrl_to_test,tmp_S];
                C_ctrl_to_test = [C_ctrl_to_test,tmp_C];
            end
        end
        p_val_uncorrected_bands_S(j,band) = ranksum(S_ctrl_to_test(:),S_exp_to_test(:));
        p_val_uncorrected_bands_C(j,band) = ranksum(C_ctrl_to_test(:),C_exp_to_test(:));
        to_bar_S{j}(1,band) = mean(S_ctrl_to_test);
        to_bar_S{j}(2,band) = mean(S_exp_to_test);
        to_bar_C{j}(1,band) = mean(C_ctrl_to_test);
        to_bar_C{j}(2,band) = mean(C_exp_to_test);
        err_S{j}(1,band) = std(S_ctrl_to_test)/sqrt(length(S_ctrl_to_test));
        err_S{j}(2,band) = std(S_exp_to_test)/sqrt(length(S_exp_to_test));
        err_C{j}(1,band) = std(C_ctrl_to_test)/sqrt(length(C_ctrl_to_test));
        err_C{j}(2,band) = std(C_exp_to_test)/sqrt(length(C_exp_to_test));
    end
end

p_val_uncorrected_bands_C = p_val_uncorrected_bands_C/2;
for i=1:length(ROI)
    [~,~,~,p_val_corrected_bands_S(i,:)] = stat_fdr_bh(p_val_uncorrected_bands_S(i,:));
    [~,~,~,p_val_corrected_bands_C(i,:)] = stat_fdr_bh(p_val_uncorrected_bands_C(i,:));
end

figure()
for region=1:length(ROI)
    subplot(1,3,region);
    p = bar(to_bar_S{region}');xticklabels(labels);
    ngroups = length(bands);
    nbars = 2;
    % Calculating the width for each bar group
    hold all
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    x = [];
    for i = 1:nbars
        x(i,:) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x(i,:), to_bar_S{region}(i,:), err_S{region}(i,:), '.');
    end
    sigstar({[x(1,1),x(2,1)],[x(1,2),x(2,2)],[x(1,3),x(2,3)],[x(1,4),x(2,4)],[x(1,5),x(2,5)],[x(1,6),x(2,6)]},p_val_corrected_bands_S(region,:)');
    title(ROI{region});ylim([-0.7,0.7]);
end
subplot(1,3,1);ylabel('Baseline normalized Spectrum');

figure();
for region=1:length(ROI)
    subplot(1,3,region);
    p = bar(to_bar_C{region}');xticklabels(labels);
    ngroups = length(bands);
    nbars = 2;
    % Calculating the width for each bar group
    hold all
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    x = [];
    for i = 1:nbars
        x(i,:) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x(i,:), to_bar_C{region}(i,:), err_C{region}(i,:), '.');
    end
    sigstar({[x(1,1),x(2,1)],[x(1,2),x(2,2)],[x(1,3),x(2,3)],[x(1,4),x(2,4)],[x(1,5),x(2,5)],[x(1,6),x(2,6)]},p_val_corrected_bands_C(region,:)');
    title(titles_c{region}(5:end));ylim([-0.8,0.8]);
end
subplot(1,3,1);ylabel('Baseline normalized Coherency');


