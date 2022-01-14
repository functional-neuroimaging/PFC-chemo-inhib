% main routine for the multielectrode LFP to calculate the spectrum and
% coherency

clear
clc
close all

%% data assignment
mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a',...
    'fr200209a','fr200327a','fr200330a','fr200331a','fr200331b'}; % list of the animals 
types = {'exp','exp','exp','exp','exp','ctrl','ctrl','ctrl','ctrl'}; % corresponding type of each animal
base_lines = {165058,162532,165221,155805,164829,165959,153031,145548,204419}; %name of the file corresponding to the end of the baseline
cnos = {172836,170719,172209,162757,172140,170359,153437,145815,204927}; % %name of the file corresponding to the beginning of CNO
n_cnos = [30,28,40,40,40,40,40,40,40]; % length of cno to analyze
n_baseline = 3; % length of baseline to analyze
fs = 1000; % sampling frequency
broken_PFC = {[31],[17,19,22,25,26,30],[],[],[],[],[],[],[]}; % list of broken channels of PFC
broken_RS = {[18,19,22,23,25,27,28,31],[13,15],[],[],[],[],[],[],[]}; % list of broken channels of Rs
broken_TH = {[],[],[4,7,15],[2,4,7,15],[4,5,7,15],[2,4,5,7,15],[],[],[]}; % list of broken channels of Th
data_folder = ''; % put the path where the preprocessed data for all the animals are there
%% load data
for MN=1:length(mouse_names)
    mouse_name = mouse_names{MN};
    path = strcat(data_folder,'\',mouse_name,'\LFP');
    type = types{MN};
    base_line_end = base_lines{MN};
    CNO_start = cnos{MN};
    n_CNO_traces = n_cnos(MN);
    brokenCh.PFC = broken_PFC{MN};
    brokenCh.Rs = broken_RS{MN};
    brokenCh.MDTh = broken_TH{MN};
    ROI = {'PFC', 'Rs', 'MDTh'};
    phase = {'baseline','CNO'};
    
    base_line_LFP = [];
    for i=1:n_baseline
        base_line_time = base_line_end - (i-1)*100;
        file_name = strcat(path,'\',mouse_name,'_',type,'_baseline_LFP_',num2str(base_line_time),'.mat');
        base_line = load(file_name);
        tmp = fieldnames(base_line);
        base_line_LFP = [base_line_LFP;base_line.(tmp{1}).dataNo50Hz];
    end
    
    clear base_line
    if size(base_line_LFP,2) == 80
        orderCh = {0, [1:32], [33:64] [65:80]}; % channel ordering
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
    
    

    cno_time = CNO_start;
    for trace=1:n_CNO_traces
        if rem(floor(cno_time/1000),10) == 6; cno_time = cno_time + 4000; end
        file_name = strcat(path,'\',mouse_name,'_',type,'_CNO_LFP_',num2str(cno_time),'.mat');
        cno = load(file_name);
        tmp = fieldnames(cno);
        cno_LFP = cno.(tmp{1}).dataNo50Hz;
        clear cno
        cno_time = cno_time + 100;
        for i=2:length(orderCh)
            LFP.CNO.(ROI{i-1}) = cno_LFP(:,orderCh{i});
        end
        clear cno_LFP
        
        for j = 1:size(ROI,2)
            LFP.CNO.(ROI{j})(:,brokenCh.(ROI{j})) = [];
        end
        
        %% Spectrogram MI
        for i=1:length(ROI)
            for j=1:length(phase)
                LFP_to_Spectrogram = mean(LFP.(phase{j}).(ROI{i}),2);
                [S,f_s,t_s] = pspectrum(LFP_to_Spectrogram, fs, 'spectrogram','FrequencyLimits',[0 100]);
                if j==1
                    Avg_factor = mean(S,2);
                else
                    S_normalized.(ROI{i}).(phase{j}) = (S - Avg_factor)./(S + Avg_factor);
                end
            end
        end
        clear Avg_factor S LFP_to_Spectrogram
        
        %% Compute Coherency MI
        win = 2000; % hyper parameters for the coherency
        noverlap = 1000;
        cxy = cell(length(phase),length(ROI));
        for i=1:length(phase)
            count = 1;
            for j=1:length(ROI)
                for k=j+1:length(ROI)
                    x = mean(LFP.(phase{i}).(ROI{j}),2);
                    y = mean(LFP.(phase{i}).(ROI{k}),2);
                    [cxy{i,count},f] = mscohere(x,y,win,noverlap,0:0.1:100,fs);
                    count = count + 1;
                end
            end
        end
        
        C_normalized = cell(length(phase)-1,3);
        titles_C = C_normalized;
        subtitles = {'PFC-Rs','PFC-MDTh','MDTh-Rs'};
        for i=1:size(cxy,2)
            for j=2:size(cxy,1)
                C_normalized{j-1,i} = (cxy{j,i} - cxy{1,i})./(cxy{j,i} + cxy{1,i});
                titles_C{j-1,i} = strcat(phase{j},'_',subtitles{i});
            end
        end
        %% save results
        save_dir = strcat(pwd,'\Results\');
        save_dir = strcat(save_dir,mouse_name(1:9),'\');
        mkdir(save_dir);
        save_file = strcat(save_dir,mouse_name(1:9),'_',num2str(trace),'_results.mat');
        save(save_file,'C_normalized','cxy','f','titles_C','S_normalized','t_s','f_s');
    end
end