% routine to process the Kir experiment, must be repeated for each animal
clear
clc 
close all

%% load data
[mouseData, mouseName, mousePath] = io_importProcessedData; % load the preprocessed Kir data

for i=1:size(mouseName,2)
    expGroup.(['d',num2str(i)]) = m_chExtractor(mouseData.(['d',num2str(i)]), 1:16);
    ctrlGroup.(['d',num2str(i)]) = m_chExtractor(mouseData.(['d',num2str(i)]), 17:32);
end

ctrlName = mouseName;
expName = mouseName;
nCh = size(ctrlGroup.d1.raw.data,2);

LFP_all_ctrl = [];
MUA_all_ctrl = [];
LFP_all_exp = [];
MUA_all_exp = [];
for i=1:size(mouseName,2)
    LFP_all_ctrl = [LFP_all_ctrl;ctrlGroup.(['d',num2str(i)]).LFP.data];
    MUA_all_ctrl = [MUA_all_ctrl;ctrlGroup.(['d',num2str(i)]).MUA.data];
    LFP_all_exp = [LFP_all_exp;expGroup.(['d',num2str(i)]).LFP.data];
    MUA_all_exp = [MUA_all_exp;expGroup.(['d',num2str(i)]).MUA.data];
end

%% Extracting Spikes 
fsMUA = ctrlGroup.d1.MUA.fs;
MUA_all_ctrl = MUA_all_ctrl(1:20*60*fsMUA,:);
MUA_all_exp = MUA_all_exp(1:20*60*fsMUA,:);
[spkTimes_ctrl, spkIndex_ctrl, ~, ~] = m_spikeDetector(MUA_all_ctrl, fsMUA);
[spkTimes_exp, spkIndex_exp, ~, ~] = m_spikeDetector(MUA_all_exp, fsMUA);
%% Spike Counting
window_size = fsMUA*60; % bin every one minute
Spike_index_ctrl = zeros(size(MUA_all_ctrl));
Spike_count_ctrl = zeros(ceil(length(MUA_all_ctrl)/window_size),size(Spike_index_ctrl,2));
Spike_index_exp = zeros(size(MUA_all_exp));
Spike_count_exp = zeros(ceil(length(MUA_all_exp)/window_size),size(Spike_index_exp,2));
for i=1:size(Spike_index_ctrl,2)
    Spike_index_ctrl(spkIndex_ctrl{i},i) = 1;
    tmp_index = Spike_index_ctrl(:,i);
    tmp_count = reshape(tmp_index,window_size,[]);
    Spike_count_ctrl(:,i) = sum(tmp_count);
    Spike_index_exp(spkIndex_exp{i},i) = 1;
    tmp_index = Spike_index_exp(:,i);
    tmp_count = reshape(tmp_index,window_size,[]);
    Spike_count_exp(:,i) = sum(tmp_count);
end
Spike_rate_exp_ch_avg = mean(Spike_count_exp,2)./60;
Spike_rate_ctrl_ch_avg = mean(Spike_count_ctrl,2)./60;
Spike_rate_change = Spike_rate_exp_ch_avg./Spike_rate_ctrl_ch_avg;
figure();plot(Spike_rate_change,'LineWidth',2);yline(1,'--r');xticks(1:2:20);ylim([0.3,1.2]);
ylabel('Control Normalized Change of Firing Rate','Fontsize',12);xlabel('Time (minutes)','FontSize',12);

save_dir = strcat(pwd,'\results'); 
mkdir(save_dir);
save_file = strcat(save_dir,'\',mouseName{1}(1:13),'_results.mat');
save(save_file,'Spike_rate_exp_ch_avg','Spike_rate_ctrl_ch_avg','Spike_rate_change','Spike_count_exp','Spike_count_ctrl');

