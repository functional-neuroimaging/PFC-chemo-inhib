%% description 
% routine to load the preprocessed data and calculate the ISI and spike
% phase relations, must be repeated for each animal


clear 
close all
clc

%% load data
[mouseData, mouseName, mousePath] = io_importProcessedData; % load required data , needs the preprocessed data
%% Spike Detection
% Concataneting MUAs 
MUA_all = [];
LFP_all = [];
analog_MUA_all = [];
time_baseline = 0;
time_CNO = 0;
fsMUA = mouseData.d1.MUA.fs;
fsLFP = mouseData.d1.LFP.fs;
for i=1:length(mouseName)
    MUA_all = [MUA_all;mouseData.(['d',num2str(i)]).MUA.data];
    LFP_all = [LFP_all;mouseData.(['d',num2str(i)]).LFP.data];
    analog_MUA_all = [analog_MUA_all;mouseData.(['d',num2str(i)]).MUA.analg_MUA];
    if contains(mouseName{i},'baseline')
       time_baseline = time_baseline + length(mouseData.(['d',num2str(i)]).MUA.data)/fsMUA;
    elseif contains(mouseName{i},'cno')
       time_CNO = time_CNO + length(mouseData.(['d',num2str(i)]).MUA.data)/fsMUA;  
    end
end
[spkTimes, ~, ~, ~] = m_spikeDetector(MUA_all, fsMUA);
clear mouseData
%% Computing ISI and phase
num_ch = length(spkTimes);
length_baseline = 15; % the duration for the baseline period in minutesn
length_transition = 15; % the duration for the transition period in minutes
length_active = 40; % the duration for the active period in minutes
end_time_baseline = time_baseline;
end_time_transition = time_baseline+(length_transition*60);
end_time_active = end_time_transition + (length_active*60);
time_LFP = linspace(0,length(LFP_all)/fsLFP,length(LFP_all));
%[b,a] = butter(3,1./(fsLFP/2),'low'); % for the infraslow band
[b,a] = butter(3,[1,4]./(fsLFP/2),'bandpass');% for the delta band
LFP_all_Low = filtfilt(b,a,LFP_all);
LFP_all_H = hilbert(LFP_all_Low);
phase_all_H = angle(LFP_all_H);

for ch=1:num_ch
    tmp = spkTimes{ch};
    phase_tmp = phase_all_H(:,ch);
    phase.baseline{ch} = phase_tmp(knnsearch(time_LFP',tmp(tmp<=end_time_baseline)));
    phase.transition{ch} = phase_tmp(knnsearch(time_LFP',tmp(tmp > end_time_baseline & tmp <= end_time_transition)));
    phase.active{ch} = phase_tmp(knnsearch(time_LFP',tmp(tmp > end_time_transition & tmp <= end_time_active)));
    spkrate.baseline{ch} = length(tmp(tmp<=end_time_baseline))./(length_baseline*60);
    spkrate.transition{ch} = length(tmp(tmp > end_time_baseline & tmp <= end_time_transition))./(length_transition*60);
    spkrate.active{ch} = length(tmp(tmp > end_time_transition & tmp <= end_time_active))./(length_active*60);
    ISI.baseline{ch} = diff(tmp(tmp<=end_time_baseline));
    ISI.transition{ch} = diff(tmp(tmp > end_time_baseline & tmp <= end_time_transition));
    ISI.active{ch} = diff(tmp(tmp > end_time_transition & tmp <= end_time_active));
end

%% save results
save_dir = strcat(pwd,'\results\MUA_statistics\');
mkdir(save_dir);
save_file = strcat(save_dir,mouseName{1}(1:13),'_results.mat');
save(save_file,'ISI','fsMUA','fsLFP','spkrate','phase');