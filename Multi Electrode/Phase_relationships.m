% routine to calculate the phase relationships between regions
clear
clc
close all
%% data assignment
mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a',...
    'fr200209a','fr200327a','fr200330a','fr200331a','fr200331b'}; % list of the animals
types = {'exp','exp','exp','exp','exp','ctrl','ctrl','ctrl','ctrl'}; % corresponding type of each animal
base_lines = {165058,162532,165221,155805,164829,165959,153031,145548,204419}; %name of the file corresponding to the end of the baseline
cnos = {174836,172519,175009,170057,175040,173459,160437,152815,211927}; % %name of the file corresponding to the beginning of CNO to analyze
broken_PFC = {[31],[17,19,22,25,26,30],[],[],[],[],[],[],[]}; % list of broken channels of PFC
broken_RS = {[18,19,22,23,25,27,28,31],[13,15],[],[],[],[],[],[],[]}; % list of broken channels of Rs
broken_TH = {[],[],[4,7,15],[2,4,7,15],[4,5,7,15],[2,4,5,7,15],[],[],[]}; % list of broken channels of Th
data_folder = ''; % put the path where the preprocessed data for all the animals are there
n_CNO_traces = 10; % the traces of cno you want yo use
n_baseline_traces = 3; % the traces of baseline you want yo use
fs = 1000; % sampling frequency of LFP
fc1 = 1;fc2 = 4; % frequency band desired
%% loading
for MN=1:length(mouse_names)
    mouse_name = mouse_names{MN};
    path = strcat(data_folder,'\',mouse_name,'\LFP');
    type = types{MN};
    base_line_end = base_lines{MN};
    CNO_start = cnos{MN};
    
    brokenCh.PFC = broken_PFC{MN};
    brokenCh.Rs = broken_RS{MN};
    brokenCh.MDTh = broken_TH{MN};
    
    ROI = {'PFC', 'Rs', 'MDTh'};
    phase = {'baseline','CNO'};
    %% load baseline LFP
    base_line_LFP = [];
    for i=1:n_baseline_traces
        base_line_time = base_line_end - (i-1)*100;
        file_name = strcat(path,'\',mouse_name,'_',type,'_baseline_LFP_',num2str(base_line_time),'.mat');
        base_line = load(file_name);
        tmp = fieldnames(base_line);
        base_line_LFP = [base_line_LFP;base_line.(tmp{1}).dataNo50Hz];
    end
    
    clear base_line
    if size(base_line_LFP,2) == 80
        orderCh = {0, [1:32], [33:64] [65:80]};
    else
        orderCh = {0, [1:32], [33:48] [49:64]};
    end
    
    for i=2:length(orderCh)
        LFP.baseline.(ROI{i-1}) = base_line_LFP(:,orderCh{i});
    end
    clear base_line_LFP
    
    for j = 1:size(ROI,2)
        LFP.baseline.(ROI{j})(:,brokenCh.(ROI{j})) = [];
    end
    
    %% load cno LFP
    cno_LFP = [];
    cno_time = CNO_start;
    for trace=1:n_CNO_traces
        if rem(floor(cno_time/1000),10) == 6; cno_time = cno_time + 4000; end
        file_name = strcat(path,'\',mouse_name,'_',type,'_CNO_LFP_',num2str(cno_time),'.mat');
        cno = load(file_name);
        tmp = fieldnames(cno);
        cno_LFP = [cno_LFP;cno.(tmp{1}).dataNo50Hz];
        clear cno
        cno_time = cno_time + 100;
    end
    
    for i=2:length(orderCh)
        LFP.CNO.(ROI{i-1}) = cno_LFP(:,orderCh{i});
    end
    clear cno_LFP
    
    for j = 1:size(ROI,2)
        LFP.CNO.(ROI{j})(:,brokenCh.(ROI{j})) = [];
    end
    
    [b1,a1] = butter(3,[fc1,fc2]./(fs/2),'bandpass');
    phase_delta_diff = cell(2,3);
    PLV = phase_delta_diff;
    for i=1:length(phase)
        for j=1:length(ROI)
            LFP_L = filtfilt(b1,a1,LFP.(phase{i}).(ROI{j}));
            anal_sig = hilbert(LFP_L);
            phase_delta.(phase{i}).(ROI{j}) = angle(anal_sig);
        end
        count = 1;
        for j=1:length(ROI)
            tmp_1 = phase_delta.(phase{i}).(ROI{j});
            for k=j+1:length(ROI)
                tmp_2 = phase_delta.(phase{i}).(ROI{k});
                count_ch = 1;
                for n_ch1=1:size(tmp_1,2)
                    for n_ch2=1:size(tmp_2,2)
                        phase_delta_diff{i,count}(:,count_ch) = tmp_1(:,n_ch1)- tmp_2(:,n_ch2);
                        ph = exp(1i*(tmp_1(:,n_ch1)- tmp_2(:,n_ch2)));
                        PLV{i,count}(:,count_ch) = abs(sum(ph))/length(ph);
                        count_ch = count_ch + 1;
                    end
                end
                count = count + 1;
            end
        end
    end
    
    save_dir = strcat(pwd,'\Results_angle\');
    mkdir(save_dir);
    save_file = strcat(save_dir,mouse_name(1:9),'_results.mat');
    save(save_file,'PLV');
end