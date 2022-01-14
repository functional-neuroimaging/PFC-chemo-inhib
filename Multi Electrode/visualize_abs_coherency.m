% routine to visualize the raw coherency, requires the results from the
% main routine, related to suuplementary figure 9

clear 
clc
close all

%% predefining
Mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a','fr200209a'...
    ,'fr200327a','fr200330a','fr200331a','fr200331b'}; % list of animals
ROI = {'PFC', 'Rs', 'MDTh'}; % recorded regions
n_exp = 5; % the first 5 are exp data
n_regions = length(ROI);
traces = {21:30,19:28,31:40,31:40,31:40,31:40,31:40,31:40,31:40}; % traces of results wants to be used
results_folder = ''; % the folder where the results from the main routine are saved

C_all_exp_base = cell(1,n_regions);
C_all_exp_cno = cell(1,n_regions);
C_all_ctrl_base = cell(1,n_regions);
C_all_ctrl_cno = cell(1,n_regions);

%% load and visualize
for reg=1:n_regions
    C_all_exp_base{reg} = [];
    C_all_ctrl_base{reg} = [];
    C_all_exp_cno{reg} = [];
    C_all_ctrl_cno{reg} = [];
    for MN=1:length(Mouse_names)
        mouse_name = Mouse_names{MN};
        path = strcat(results_folder,'\',mouse_name);
        n_traces = traces{MN};
        C_temp = [];
        for trace=n_traces(1):n_traces(end)
            file_name = strcat(path,'\',mouse_name,'_',num2str(trace),'_results','.mat');
            Data = load(file_name);
            C_temp = [C_temp;Data.cxy{2,reg}];
        end
        if MN<=n_exp
            C_all_exp_cno{reg} = [C_all_exp_cno{reg};C_temp];
            C_all_exp_base{reg} = [C_all_exp_base{reg};Data.cxy{1,reg}];
        else
            C_all_ctrl_cno{reg} = [C_all_ctrl_cno{reg};C_temp];
            C_all_ctrl_base{reg} = [C_all_ctrl_base{reg};Data.cxy{1,reg}];
        end
        fc = Data.f;
    end
    subplot(2,n_regions,reg);title(Data.titles_C{reg}(5:end));
    to_plot_base = mean(C_all_exp_base{reg});
    err_base = std(C_all_exp_base{reg})/sqrt(size(C_all_exp_base{reg},1));
    plot_shadedErrorBar(fc(1:991),to_plot_base(1:991),err_base(1:991),'lineProps','-b');
    hold all
    to_plot_cno = mean(C_all_exp_cno{reg});
    err_cno = std(C_all_exp_cno{reg})/sqrt(size(C_all_exp_cno{reg},1));
    plot_shadedErrorBar(fc(1:991),to_plot_cno(1:991),err_cno(1:991),'lineProps','-r');
    subplot(2,3,reg+3);
    to_plot_base = mean(C_all_ctrl_base{reg});
    err_base = std(C_all_ctrl_base{reg})/sqrt(size(C_all_ctrl_base{reg},1));
    plot_shadedErrorBar(fc(1:991),to_plot_base(1:991),err_base(1:991),'lineProps','-b');
    hold all
    to_plot_cno = mean(C_all_ctrl_cno{reg});
    err_cno = std(C_all_ctrl_cno{reg})/sqrt(size(C_all_ctrl_cno{reg},1));
    plot_shadedErrorBar(fc(1:991),to_plot_cno(1:991),err_cno(1:991),'lineProps','-r');
end
subplot(2,3,1);ylabel('abs of coherency');
subplot(2,3,4);ylabel('abs of coherency');
subplot(2,3,4);xlabel('frequency (Hz)');
subplot(2,3,5);xlabel('frequency (Hz)');
subplot(2,3,6);xlabel('frequency (Hz)');