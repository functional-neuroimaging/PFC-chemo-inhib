% routine to calculate the gamma enevelope based on Nir paper
clear
clc
close all
%% data assignment
mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a',...
    'fr200209a','fr200327a','fr200330a','fr200331a','fr200331b'}; % list of the animals
types = {'exp','exp','exp','exp','exp','ctrl','ctrl','ctrl','ctrl'}; % corresponding type of each animal
base_lines = {165058,162532,165221,155805,164829,165959,153031,145548,204419}; %name of the file corresponding to the end of the baseline
cnos = {174836,172519,175009,170057,175040,173459,160437,152815,211927}; %name of the file corresponding to the beginning of CNO to analyzed
broken_PFC = {[31],[17,19,22,25,26,30],[],[],[],[],[],[],[]}; % list of broken channels of PFC
broken_RS = {[18,19,22,23,25,27,28,31],[13,15],[],[],[],[],[],[],[]}; % list of broken channels of Rs
broken_TH = {[],[],[4,7,15],[2,4,7,15],[4,5,7,15],[2,4,5,7,15],[],[],[]}; % list of broken channels of Th
data_folder = ''; % put the path of the required data
n_CNO_traces = 10; % number of traces of cno you want yo use
n_baseline_traces = 3; % number of traces of baseline you want yo use
fs = 1000; % sampling frequency
fc1 = 30;fc2 = 70; % gamma frequency
ROI = {'PFC', 'Rs', 'MDTh'}; % regions recorded
phase = {'baseline','CNO'}; % different phases of recording
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
        
        for i=2:length(orderCh)
            LFP.CNO.(ROI{i-1}) = cno_LFP(:,orderCh{i});
        end
        clear cno_LFP
        
        for j = 1:size(ROI,2)
            LFP.CNO.(ROI{j})(:,brokenCh.(ROI{j})) = [];
        end
        
        %% calculating gamma envelope
        envelope_gamma = [];
        for ph=1:length(phase)
            for reg=1:length(ROI)
                LFP_to_inves = mean(LFP.(phase{ph}).(ROI{reg}),2);
                [S,f_spectro,t_spec] = pspectrum(LFP_to_inves, fs, 'spectrogram','FrequencyLimits',[0 80],...
                    'timeResolution',2,'OverlapPercent',75,'Leakage',0.85);
                S_gamma = log10(S(f_spectro>30 & f_spectro<=70,:));
                envl = (sum(S_gamma)) - mean((sum(S_gamma)));
                envelope_gamma.(phase{ph}).(ROI{reg}) = envl;
                fs_envelope = 1./(t_spec(2) - t_spec(1));
                [p_envelope.(phase{ph}).(ROI{reg}),fspec_envelope] = pspectrum(envelope_gamma.(phase{ph}).(ROI{reg}),fs_envelope,'frequencyresolution',0.05,'frequencylimit',[0.001,1]);
            end
        end
        
        cxy_envelope = cell(length(phase),length(ROI));
        for i=1:length(phase)
            count = 1;
            for j=1:length(ROI)
                for k=j+1:length(ROI)
                    x = envelope_gamma.(phase{i}).(ROI{j});
                    y = envelope_gamma.(phase{i}).(ROI{k});
                    [cxy_envelope{i,count},fc_envelope] = mscohere(x,y,[],[],0.001:0.004:1,fs_envelope);
                    count = count + 1;
                end
            end
        end
        
        C_normalized_envelope = cell(length(phase)-1,3);
        for i=1:size(cxy_envelope,2)
            for j=2:size(cxy_envelope,1)
                C_normalized_envelope{j-1,i} = (cxy_envelope{j,i} - cxy_envelope{1,i})./(cxy_envelope{j,i} + cxy_envelope{1,i});
            end
        end
        
        save_dir = strcat(pwd,'\Results_gamma_envelope\');
        save_dir = strcat(save_dir,mouse_name(1:9),'\');
        mkdir(save_dir);
        save_file = strcat(save_dir,mouse_name(1:9),'_',num2str(trace),'_results.mat');
        save(save_file,'cxy_envelope','C_normalized_envelope',...
            'fc_envelope','p_envelope','fspec_envelope');
    end
end