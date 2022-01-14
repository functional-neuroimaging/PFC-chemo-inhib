% routine to visualize the raw coherency for the gamma envelope, requires the results from the
% Gamma_envelope routine, related to suuplementary figure 9


clear all
clc 
close all

%% predefining
Mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a','fr200209a'...
    ,'fr200327a','fr200330a','fr200331a','fr200331b'}; % list of animals
ROI = {'PFC-RS','PFC-MDTh','RS-MDTh'}; % name of the pair of the regions
results_folder = '';% put the folder where the results are
n_exp = 5; % the first 5 are the exp data

C_all_exp_base = cell(1,length(ROI));
C_all_exp_cno = cell(1,length(ROI));
C_all_ctrl_base = cell(1,length(ROI));
C_all_ctrl_cno = cell(1,length(ROI));
%% load results and visualize
for reg=1:length(ROI)
    C_all_exp_base{reg} = [];
    C_all_ctrl_base{reg} = [];
    C_all_exp_cno{reg} = [];
    C_all_ctrl_cno{reg} = [];
    for MN=1:length(Mouse_names)
        mouse_name = Mouse_names{MN};
        file_name = strcat(results_folder,'\',mouse_name,'_results','.mat');
        Data = load(file_name);
        if MN<=n_exp
            C_all_exp_cno{reg} = [C_all_exp_cno{reg};Data.cxy_envelope{2,reg}];
            C_all_exp_base{reg} = [C_all_exp_base{reg};Data.cxy_envelope{1,reg}];
        else
            C_all_ctrl_cno{reg} = [C_all_ctrl_cno{reg};Data.cxy_envelope{2,reg}];
            C_all_ctrl_base{reg} = [C_all_ctrl_base{reg};Data.cxy_envelope{1,reg}];
        end
        fc = Data.fc_envelope;
    end
    subplot(2,3,reg);title(Data.titles_C{reg}(5:end));
    to_plot_base = mean(C_all_exp_base{reg});
    err_base = std(C_all_exp_base{reg})/sqrt(5);
    plot_shadedErrorBar(fc(fc<0.5),to_plot_base(fc<0.5),err_base(fc<0.5),'lineProps','-b');
    ylim([0,1]);
    hold all
    to_plot_cno = mean(C_all_exp_cno{reg});
    err_cno = std(C_all_exp_cno{reg})/sqrt(5);
    plot_shadedErrorBar(fc(fc<0.5),to_plot_cno(fc<0.5),err_cno(fc<0.5),'lineProps','-r');
    ylim([0,1]);
    subplot(2,3,reg+3);
    to_plot_base = mean(C_all_ctrl_base{reg});
    err_base = std(C_all_ctrl_base{reg})/sqrt(5);
    plot_shadedErrorBar(fc(fc<0.5),to_plot_base(fc<0.5),err_base(fc<0.5),'lineProps','-b');
    ylim([0,1]);
    hold all
    to_plot_cno = mean(C_all_ctrl_cno{reg});
    err_cno = std(C_all_ctrl_cno{reg})/sqrt(5);
    plot_shadedErrorBar(fc(fc<0.5),to_plot_cno(fc<0.5),err_cno(fc<0.5),'lineProps','-r');
    ylim([0,1]);
end
subplot(2,3,1);ylabel('abs of coherency');
subplot(2,3,4);ylabel('abs of coherency');
subplot(2,3,4);xlabel('frequency (Hz)');
subplot(2,3,5);xlabel('frequency (Hz)');
subplot(2,3,6);xlabel('frequency (Hz)');