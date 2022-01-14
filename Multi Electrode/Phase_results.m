% routine to check the results for phase, requires the
% results from the phase_relationships routine, related to figure 5

clear
clc
close all

%% load results
[file,path] = uigetfile('*.mat','Select One or More results data', ...
    'MultiSelect', 'on'); % load the results from the phase_relationships
Data = cell(1,length(file));
mousename = cell(1,length(file));
for i=1:length(file)
    Data{i} = load([path,file{i}]);
    mousename{i} = file{i}(1:end-4);
end

ROI = {'PFC', 'Rs', 'MDTh'}; % regions recorded
n_exp = 5; % the first 5 are exp

PLV_all_exp_between = cell(2,length(ROI));
PLV_all_ctrl_between = cell(2,length(ROI));

for i=1:length(Data)
    for j=1:length(ROI)
        if i<=n_exp
            PLV_all_exp_between{1,j} = [PLV_all_exp_between{1,j},Data{i}.PLV{1,j}];
            PLV_all_exp_between{2,j} = [PLV_all_exp_between{2,j},Data{i}.PLV{2,j}];
        else
            PLV_all_ctrl_between{1,j} = [PLV_all_ctrl_between{1,j},Data{i}.PLV{1,j}];
            PLV_all_ctrl_between{2,j} = [PLV_all_ctrl_between{2,j},Data{i}.PLV{2,j}];
        end
    end
end

%% Modulation index
MI_exp = cell(1,length(ROI));
MI_ctrl = cell(1,length(ROI));
to_bar_MI = zeros(2,length(ROI));
se_MI = zeros(2,length(ROI));
p_val_MI = zeros(1,length(ROI));
for i=1:length(ROI)
    MI_exp{i} = (PLV_all_exp_between{2,i} - PLV_all_exp_between{1,i})./(PLV_all_exp_between{2,i} + PLV_all_exp_between{1,i});
    MI_ctrl{i} = (PLV_all_ctrl_between{2,i} - PLV_all_ctrl_between{1,i})./(PLV_all_ctrl_between{2,i} + PLV_all_ctrl_between{1,i});
    to_bar_MI(:,i) = [mean(MI_ctrl{i}),mean(MI_exp{i})];
    se_MI(:,i) = [std(MI_ctrl{i}),std(MI_exp{i})]./[sqrt(length(MI_ctrl{i})),sqrt(length(MI_exp{i}))];
    p_val_MI(i) = ranksum(MI_ctrl{i},MI_exp{i});
end
corrected_p_val_MI = bonf_holm(p_val_MI);

x = [];
ngroups = length(ROI);
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x(:,i) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
end
groups ={[x(1,1),x(1,2)],[x(2,1),x(2,2)],[x(3,1),x(3,2)]};
figure()
bar(to_bar_MI');ylabel('PLV Modulation Index')
xticklabels({'PFC-RS','PFC-MDTh','Rs-MDTh'});legend({'Control','Experimental'},'Location','northeastoutside');
hold on
er = errorbar(x(1,:),to_bar_MI(:,1),se_MI(:,1));
er.Color = [0 0 0];er.LineStyle = 'none';  
er = errorbar(x(2,:),to_bar_MI(:,2),se_MI(:,2));
er.Color = [0 0 0];er.LineStyle = 'none'; 
er = errorbar(x(3,:),to_bar_MI(:,3),se_MI(:,3));
er.Color = [0 0 0];er.LineStyle = 'none'; 
sigstar(groups,corrected_p_val_MI)