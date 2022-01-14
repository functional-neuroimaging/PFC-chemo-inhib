% routine to check the results for the spike rate, needs the results from
% main routine, related to supplementary figure 1

clear 
clc
close all

%% load results 
[file,path] = uigetfile('*.mat','Select One or More results data', ...
    'MultiSelect', 'on'); % load the results derived from the main routine
Data = cell(1,length(file));
mousename = cell(1,length(file));
for i=1:length(file)
    Data{i} = load([path,file{i}]);
    mousename{i} = file{i}(1:end-4);
end

%% check rate
n_animals = 4;
length_data = 20; % data length in minutes
R_exp = zeros(1,n_animals);
R_ctrl = zeros(1,n_animals);
for i=1:length(Data)
    R_exp(i) = mean(sum(Data{i}.Spike_count_exp))/(length_data*60);
    R_ctrl(i) = mean(sum(Data{i}.Spike_count_ctrl))/(length_data*60);
end
[~,p_val_all] = ttest(R_ctrl,R_exp,'Tail','right');
figure()
for i=1:n_animals
    plot(1,R_ctrl(i),'xb','LineWidth',4)
    hold all
    plot(2,R_exp(i),'xb','LineWidth',4)
end
plot([R_ctrl;R_exp],'k','LineWidth',0.5)
xticks([1,2]);xticklabels({'Ctrl','Exp'});xlim([0.8,2.2]);ylabel('Channel Averaged Firing Rate','FontSize',12);
sigstar({[1,2]},p_val_all)