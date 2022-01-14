% routine to check the results for gamma envelope coherency, requires the
% results from the BLP_gamma_envelope routine

clear
clc
close all

%% load results
Mouse_names = {'fr191209a','fr200121a','fr200207a','fr200208a','fr200209a','fr200327a'...
    ,'fr200330a','fr200331a','fr200331b'}; % name of the animals
results_folder = 'E:\Noei\PhD\Projects\Carola_Project\Coherence_proj\Results_gamma_envelope_pooling'; % put the directory where the results are saved
C_pool_one = cell(1,length(Mouse_names));
n_traces = 10; % number of traces you want to load
ROI = {'PFC', 'Rs', 'MDTh'};
for MN=1:length(Mouse_names)
    mouse_name = Mouse_names{MN};
    path = strcat(results_folder,'\',mouse_name);
    
    C_pool_one{MN} = cell(length(ROI),1);
    C_all = cell(1,length(ROI));
    
    for trace=1:n_traces
        file_name = strcat(path,'\',mouse_name,'_',num2str(trace),'_results','.mat');
        Data = load(file_name);
        
        for i=1:length(ROI)
            C_all{i} = [C_all{i};Data.C_normalized_envelope{i}];
            alpha_c = Data.C_normalized_envelope{i};
            C_pool_one{MN}{i} = [C_pool_one{MN}{i};alpha_c];
        end
    end
    
    
    C_downsampled = C_all;
    C_pool_all{MN} = C_downsampled;
    fc_all{MN} = Data.fc_envelope;
end
titles_c = {'PFC-Rs','PFC-MDTh','Rs-MDTh'};
bands = {[0,0.5],[0.5,1]};
labels = {'infraslow','slow'};


for j=1:length(ROI)
    for band=1:length(bands)
        C_exp_to_test = [];
        C_ctrl_to_test = [];
        for i=1:length(Mouse_names)
            tmp_C = C_pool_one{i}{j}(:,fc_all{MN} > bands{band}(1) & fc_all{MN} <= bands{band}(2));
            tmp_C = median(tmp_C,2);
            tmp_C = tmp_C';
            if i<=5
                C_exp_to_test = [C_exp_to_test,tmp_C];
            else
                C_ctrl_to_test = [C_ctrl_to_test,tmp_C];
            end
        end
        p_val_uncorrected_bands_C(j,band) = ranksum(C_ctrl_to_test(:),C_exp_to_test(:));
        to_bar_C{j}(1,band) = mean(C_ctrl_to_test);
        to_bar_C{j}(2,band) = mean(C_exp_to_test);
        err_C{j}(1,band) = std(C_ctrl_to_test)/sqrt(length(C_ctrl_to_test));
        err_C{j}(2,band) = std(C_exp_to_test)/sqrt(length(C_exp_to_test));
    end
end

p_val_uncorrected_bands_C = p_val_uncorrected_bands_C./2;
for i=1:length(ROI)
    [~,~,~,p_val_corrected_bands_C(i,:)] = stat_fdr_bh(p_val_uncorrected_bands_C(i,:));
end


figure();
for region=1:length(ROI)
    subplot(1,length(ROI),region);%title(ROI{region});
    p = bar(to_bar_C{region}');xticklabels(labels);
    ngroups = 2;
    nbars = 2;
    % Calculating the width for each bar group
    hold all
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    x = [];
    for i = 1:nbars
        x(i,:) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x(i,:), to_bar_C{region}(i,:), err_C{region}(i,:), '.');
    end
    sigstar({[x(1,1),x(2,1)],[x(1,2),x(2,2)]},p_val_corrected_bands_C(region,:)');
    title(titles_c{region});ylim([-0.8,0.8]);
    %ax = gca;ax.FontSize = 30;
    %legend([p(1),p(2)],'Ctrl','Exp')
end
subplot(1,length(ROI),1);ylabel('Baseline normalized Coherency');