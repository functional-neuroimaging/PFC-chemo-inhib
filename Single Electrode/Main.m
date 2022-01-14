%% description 
% routine to load the preprocessed data and calculate the spike rate and
% normalized spectrograms, must be done for each animal individually


clear 
close all
clc


%% load data
[mouseData, mouseName, mousePath] = io_importProcessedData; % load required data , needs the preprocessed data
ch = 0; % 1 for computing channel wise, 0 for averaging over channels
normalize = 1; % 0 for the case we want to use only the last 5 minutes, 1 for all the baseline
th = 1; % set 1 for th on the whole data, 0 for the th on the baseline
all = 1; % 1 if you want to save everything
%% Base line CNO extraction
% finding the 5 minute baseline streams
correct_baseline = zeros(1,length(mouseName));
correct_CNO = zeros(1,length(mouseName));
for i=1:length(mouseName)
    fs = mouseData.(['d',num2str(i)]).LFP.fs;
    correct_baseline(i) = (fs*300 == length(mouseData.(['d',num2str(i)]).LFP.data) && contains(mouseName{i},'baseline'));
    correct_CNO(i) = (fs*300 == length(mouseData.(['d',num2str(i)]).LFP.data) && contains(mouseName{i},'cno'));
end
idx_correct_baseline = find(correct_baseline == 1);
idx_correct_cno = find(correct_CNO == 1);
%% Normalized LFP Spectrogram
% computing average baseline Spectrogram
S_baseline = [];

if normalize == 1 && ch == 0
    for i=1:length(idx_correct_baseline)
        fs = mouseData.(['d',num2str(idx_correct_baseline(i))]).LFP.fs;
        Data_to_spectrogram = mean(mouseData.(['d',num2str(idx_correct_baseline(i))]).LFP.data,2);
        [S_baseline(i,:,:),~,~] = pspectrum(Data_to_spectrogram, fs, 'spectrogram','FrequencyLimits',[0 150]);
    end
    S_baseline_avg = squeeze(mean(S_baseline,1));
    Avg_factor = squeeze(mean(S_baseline_avg,2));
elseif normalize == 0
    fs = mouseData.(['d',num2str(idx_correct_baseline(end))]).LFP.fs;
    Data_to_spectrogram = mean(mouseData.(['d',num2str(idx_correct_baseline(end))]).LFP.data,2);
    [S_baseline_avg,~,~] = pspectrum(Data_to_spectrogram, fs, 'spectrogram','FrequencyLimits',[0 150]);
    Avg_factor = squeeze(mean(S_baseline_avg,2));
elseif normalize == 1 && ch == 1 
    for i=1:length(idx_correct_baseline)
        fs = mouseData.(['d',num2str(idx_correct_baseline(i))]).LFP.fs;
        Data_to_spectrogram = mouseData.(['d',num2str(idx_correct_baseline(i))]).LFP.data;
        for j=1:size(Data_to_spectrogram,2)
            [S_baseline(i,j,:,:),~,~] = pspectrum(Data_to_spectrogram(:,j), fs, 'spectrogram','FrequencyLimits',[0 150]);
            S_baseline_avg = squeeze(mean(S_baseline,1));
        end
    end
    Avg_factor = squeeze(mean(S_baseline_avg,3));
end


% computing normalized spectrograms
all_cno_LFP = [];

for i = 1:length(idx_correct_cno)
    if ch == 0
        all_cno_LFP = [all_cno_LFP;mean(mouseData.(['d',num2str(idx_correct_cno(i))]).LFP.data,2)];
    else
        all_cno_LFP = [all_cno_LFP;mouseData.(['d',num2str(idx_correct_cno(i))]).LFP.data];
    end
end

S_all = [];
for i=1:size(all_cno_LFP,2)
    [S_all(i,:,:),f,t] = pspectrum(all_cno_LFP(:,i), fs, 'spectrogram','FrequencyLimits',[0 150]);
end
S_all = squeeze(S_all);
S_all_normalized = (S_all - Avg_factor)./(S_all + Avg_factor);
t_scale = t/60;
%% Spike Detection
% Concataneting MUAs 
MUA_all = [];
for i=1:length(mouseName)
    if correct_baseline(i) || correct_CNO(i)
        MUA_all = [MUA_all;mouseData.(['d',num2str(i)]).MUA.data];
    end
end
fsMUA = mouseData.d1.MUA.fs;
if th == 1
    [spkTimes, spkIndex, ~, ~] = m_spikeDetector(MUA_all, fsMUA);
else
    MUA_baseline = [];
    for i=1:length(idx_correct_baseline)
        MUA_baseline = [MUA_baseline;mouseData.(['d',num2str(idx_correct_baseline(i))]).MUA.data];
    end
    [~,~, ~,spk_th] = m_spikeDetector(MUA_baseline, fsMUA);
    [spkTimes, spkIndex, ~, ~] = m_spikeDetector(MUA_all, fsMUA,[],spk_th);
end
%% Spike Counting
Spike_index_all = zeros(size(MUA_all));
window_size = fsMUA*60; % bin every one minute
Spike_count = zeros(ceil(length(MUA_all)/window_size),size(Spike_index_all,2));
for i=1:size(Spike_index_all,2)
    Spike_index_all(spkIndex{i},i) = 1;
    tmp_index = Spike_index_all(:,i);
    tmp_count = reshape(tmp_index,window_size,[]);
    Spike_count(:,i) = sum(tmp_count);
end
Spike_rate = Spike_count./60;
Spike_rate_ch_avg = mean(Spike_count,2)./60;
idx_baseline_end = length(idx_correct_baseline)*5;
if normalize == 1 && ch == 0
    baseline_avg_Spik_rate_ch_avg = mean(Spike_rate_ch_avg(1:idx_baseline_end));
    cno_spike_rate_percent_change_ch_avg = Spike_rate_ch_avg(idx_baseline_end+1:end)./baseline_avg_Spik_rate_ch_avg;
elseif normalize == 1 && ch == 1
    baseline_avg_Spik_rate = mean(Spike_rate(1:idx_baseline_end,:),1);
    cno_spike_rate_percent_change = Spike_rate(idx_baseline_end+1:end,:)./baseline_avg_Spik_rate;
else
    baseline_avg_Spik_rate_ch_avg = mean(Spike_rate_ch_avg(idx_baseline_end-1:idx_baseline_end));
    cno_spike_rate_percent_change_ch_avg = Spike_rate_ch_avg(idx_baseline_end+1:end)./baseline_avg_Spik_rate_ch_avg;
end
%figure();plot(cno_spike_rate_percent_change_ch_avg);yline(1,'--r');
%% Save results
%save_dir = 'E:\Noei\PhD\Projects\Carola_Project\Data\Results_Shahryar\';
save_dir = strcat(pwd,'\results\');
if th == 1 && normalize == 1 && ch == 0
    tmp_str = 'All_baseline-Th_whole\';
elseif th == 1 && normalize == 0
    tmp_str = 'L5_baseline_Th_whole\';
elseif th == 0 && normalize == 1
    tmp_str = 'All_baseline_Th_baseline\';
elseif th == 0 && normalize == 0
    tmp_str = 'L5_baseline_Th_baseline\';
elseif th == 1 && normalize == 1 && ch == 1
    tmp_str = 'All_baseline-Th_whole-chlvl\';
elseif all == 1
    tmp_str = 'All_baseline-Th_whole-time_resolved\';
end
save_dir = strcat(save_dir,tmp_str);
mkdir(save_dir)
save_file = strcat(save_dir,mouseName{1}(1:13),'_results.mat');
if ch == 0 && all == 0
    save(save_file,'cno_spike_rate_percent_change_ch_avg','S_all_normalized','t','f');
elseif ch == 1 && all == 0
    save(save_file,'cno_spike_rate_percent_change','S_all_normalized','S_all','S_baseline_avg','t','f');
elseif ch ==0 && all == 1
    save(save_file,'cno_spike_rate_percent_change_ch_avg','S_all_normalized','S_all','S_baseline','t','f');
end