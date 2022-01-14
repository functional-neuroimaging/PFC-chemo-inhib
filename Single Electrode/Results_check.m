% routine to visualize and calculate the statistiscs, requires the resulsts obtained from the main
% and the spiking properties, related to the figure 2 and figure 4

close all
clear 
clc

%% load results 
[file,path] = uigetfile('*.mat','Select One or More results data', ...
    'MultiSelect', 'on'); % load the results that are obtained from the main routine
Data = cell(1,length(file));
mousename = cell(1,length(file));
for i=1:length(file)
    Data{i} = load([path,file{i}]);
    mousename{i} = file{i}(1:end-4);
end

%% exp/ctrl rate comparison 

% seperating exp and ctrl
exp_idx = find(contains(mousename,'exp'));
rate_all_exp = zeros(55,length(exp_idx));
for i=1:length(exp_idx)
    rate_all_exp(:,i) = Data{exp_idx(i)}.cno_spike_rate_percent_change_ch_avg(1:55);
end

ctr_idx = find(contains(mousename,'ctr'));
rate_all_ctr = zeros(55,length(ctr_idx));
for i =1:length(ctr_idx)
    rate_all_ctr(:,i) = Data{ctr_idx(i)}.cno_spike_rate_percent_change_ch_avg(1:55);
end

% visaulizing overtime
for i=1:size(rate_all_exp,1)
    [p_vs_baseline(i),h_vs_baseline(i)] = signrank(rate_all_exp(i,:)',1,'Tail','left');
    [p_vs_ctrl(i),h_vs_ctrl(i)] = ranksum(rate_all_exp(i,:)',rate_all_ctr(i,:)','Tail','left');
    [p_ctrl_vs_baseline(i),h_ctrl_vs_baseline(i)] = signrank(rate_all_ctr(i,:)',1,'Tail','left');
end
h_vs_baseline_corrected = stat_fdr_bh(p_vs_baseline);
h_vs_ctrl_corrected = stat_fdr_bh(p_vs_ctrl);

figure();h1 = plot_shadedErrorBar([],mean(rate_all_exp,2),std(rate_all_exp,[],2)/sqrt(size(rate_all_exp,2)));h1.patch.FaceColor = 'b';h1.mainLine.Color = 'b';h1.mainLine.LineWidth = 2;
hold on;h2 = plot_shadedErrorBar([],mean(rate_all_ctr,2),std(rate_all_ctr,[],2)/sqrt(size(rate_all_ctr,2)));h2.patch.FaceColor = 'g';h2.mainLine.Color = 'g';h2.mainLine.LineWidth = 2;
yline(1,'--k');ylim([0,2]);plot(find(h_vs_ctrl),0,'x','Color','r','LineWidth',4);xlim([0,55])
xlabel('Time after injection (minutes)','FontSize',12);ylabel('Baseline normalized rate','FontSize',12);axis square
legend([h1.mainLine,h2.mainLine],{'Exp','Ctrl'},'Location','northeast','Box','off','FontSize',8)

% checking for autocorrelation among data
p_val_uncorrected_gaussian = zeros(1,2*size(rate_all_ctr,2));
lag_exp = zeros(1,size(rate_all_exp,2));
lag_ctr = zeros(1,size(rate_all_ctr,2));
count = 1;
for i=1:size(rate_all_exp,2)
    [~, p_val_uncorrected_gaussian(count), ~] = stat_swtest(rate_all_exp(:,i));
    count = count+1;
    [~, p_val_uncorrected_gaussian(count), ~] = stat_swtest(rate_all_ctr(:,i));
    count = count+1;
    [acf,~,bounds] = autocorr(rate_all_exp(:,i));
    lag_exp(i) = find(abs(acf) < bounds(1),1,'first');
    [acf,~,bounds] = autocorr(rate_all_ctr(:,i));
    lag_ctr(i) = find(abs(acf) < bounds(1),1,'first');
end
h_gaussian = stat_fdr_bh(p_val_uncorrected_gaussian);

% pooling based on independency
rate_decreased_exp = cell(1,size(rate_all_exp,2));
time_decreased_exp = rate_decreased_exp;
rate_decreased_ctrl = cell(1,size(rate_all_ctr,2));
time_decreased_ctrl = rate_decreased_ctrl;
for i=1:size(rate_all_exp,2) 
    rate_decreased_exp{i} = rate_all_exp(1:lag_exp(i):end,i);
    time_decreased_exp{i} = 1:lag_exp(i):55;
    rate_decreased_ctrl{i} = rate_all_ctr(1:lag_ctr(i):end,i);
    time_decreased_ctrl{i} = 1:lag_ctr(i):55;
end

pooling_idx = [1,15;15,35;35,55];
to_plot_ctrl = zeros(1,size(pooling_idx,1));
to_plot_exp = zeros(1,size(pooling_idx,1));
error_bar_ctrl = zeros(1,size(pooling_idx,1));
error_bar_exp = zeros(1,size(pooling_idx,1));
p_val_uncorrected_pool = zeros(1,size(pooling_idx,1));
for i=1:size(pooling_idx,1)
    to_test_exp = [];
    to_test_ctrl = [];
    for j=1:size(rate_decreased_ctrl,2)
        tmp_ctrl = rate_decreased_ctrl{j};
        tmp_exp = rate_decreased_exp{j};
        time_ctrl = time_decreased_ctrl{j};
        time_exp  = time_decreased_exp{j};
        to_test_exp = [to_test_exp;tmp_exp(time_exp >= pooling_idx(i,1) & time_exp < pooling_idx(i,2)-1)];
        to_test_ctrl = [to_test_ctrl;tmp_ctrl(time_ctrl >= pooling_idx(i,1)& time_ctrl < pooling_idx(i,2)-1)];
    end
    p_val_uncorrected_pool(i) = ranksum(to_test_exp(:),to_test_ctrl(:),'tail','left');
    to_plot_ctrl(i) = mean(to_test_ctrl(:));
    to_plot_exp(i) = mean(to_test_exp(:));
    error_bar_ctrl(i) = std(to_test_ctrl(:))/sqrt(length(to_test_ctrl(:)));
    error_bar_exp(i) = std(to_test_exp(:))/sqrt(length(to_test_exp(:)));
end
p_val_corrected_pool = bonf_holm(p_val_uncorrected_pool);
figure();hndl_err_ctrl = errorbar(to_plot_ctrl,error_bar_ctrl,'g','LineWidth',2,'HandleVisibility','off');hold all;sigstar({1,2,3},p_val_corrected_pool);
hndl_err_exp = errorbar(to_plot_exp,error_bar_exp,'b','LineWidth',2);ylim([0.2,1.2]);xlim([0,4]);axis square;
xticks(1:6);xticklabels({'0-15','15-35','35-55'});xtickangle(45);
xlabel('Time after injection (minutes)','FontSize',12);ax = ancestor(gca, 'axes');yrule = ax.YAxis;yrule.FontSize = 8;
ylabel('Baseline normalized rate','FontSize',12);
legend([hndl_err_exp,hndl_err_ctrl],{'exp','ctrl'},'Location','northeast','Box','off','FontSize',8);

%% exp/ctrl spectrogram comparison
S_exp_all = [];
time_exp_minutes = [];
freq_exp_Hz = [];
for i=1:length(exp_idx)
    S_exp_all(:,:,i) = Data{exp_idx(i)}.S_all_normalized;
    time_exp_minutes(:,i) = Data{exp_idx(i)}.t./60;
    freq_exp_Hz(:,i) = Data{exp_idx(i)}.f;
end

S_ctr_all = [];
freq_ctr_Hz = [];
time_ctr_minutes = [];
for i =1:length(ctr_idx)
    S_ctr_all(:,:,i) = Data{ctr_idx(i)}.S_all_normalized;
    time_ctr_minutes(:,i) = Data{ctr_idx(i)}.t./60;
    freq_ctr_Hz(:,i) = Data{ctr_idx(i)}.f;
end

% pooling based on autocorr
delta_time_exp = cell(1,size(rate_all_exp,2));
delta_time_ctrl = cell(1,size(rate_all_ctr,2));
count = 1;
p_val_uncorrected_gaussian_freq = [];
S_exp_downsampled = delta_time_exp;
S_ctrl_downsampled = delta_time_ctrl;
time_ctrl_dowmsampled = delta_time_ctrl;
time_exp_dowmsampled = delta_time_exp;
for i=1:length(delta_time_exp)
    S_tmp_exp = squeeze(S_exp_all(:,:,i));
    S_tmp_ctrl = squeeze(S_ctr_all(:,:,i));
    lag_exp = zeros(1,size(S_tmp_exp,1));
    lag_ctrl = zeros(1,size(S_tmp_ctrl,1));
    for j=1:size(S_tmp_exp,1)
        [acf,~,bounds] = autocorr(S_tmp_exp(j,:),200);
        [~,p_val_uncorrected_gaussian_freq(count)] = stat_swtest(S_tmp_exp(j,:));
        count = count + 1;
        lag_exp(j) = find(abs(acf) < bounds(1),1,'first');
        [acf,~,bounds] = autocorr(S_tmp_ctrl(j,:),200);
        [~,p_val_uncorrected_gaussian_freq(count)] = stat_swtest(S_tmp_ctrl(j,:));
        count = count + 1;
        lag_ctrl(j) = find(abs(acf) < bounds(1),1,'first');
    end
    delta_time_exp{i} = floor(mean(lag_exp));
    delta_time_ctrl{i} = floor(mean(lag_ctrl));
    S_exp_downsampled{i} = S_tmp_exp(:,1:delta_time_exp{i}:end);
    S_ctrl_downsampled{i} = S_tmp_ctrl(:,1:delta_time_ctrl{i}:end);
    time_exp_dowmsampled{i} = time_exp_minutes(1:delta_time_exp{i}:end,i);
    time_ctrl_dowmsampled{i} = time_ctr_minutes(1:delta_time_ctrl{i}:end,i);
end
h_gaussian_freq = stat_fdr_bh(p_val_uncorrected_gaussian_freq);

% over bands
bands = {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70],[70,100]};
to_bar = zeros(2,length(bands));
err = to_bar;
labels = {'infraslow','Delta','Theta','Alpha','Beta','Low Gamma','High Gamma'};
p_val_uncorrected_bands = zeros(1,length(bands));
for i=1:length(bands)
    freq_ctr_idx = freq_ctr_Hz(:,1) >= bands{i}(1) & freq_ctr_Hz(:,1) < bands{i}(2);
    freq_exp_idx = freq_ctr_Hz(:,1) >= bands{i}(1) & freq_ctr_Hz(:,1) < bands{i}(2);
    S_ctr_to_test = [];
    S_exp_to_test = [];
    for k=1:size(S_ctr_all,3)
        tmp_ctrl = S_ctrl_downsampled{k};
        time_ctr_idx = time_ctrl_dowmsampled{k} >= 15 & time_ctrl_dowmsampled{k} <=55;
        tmp_ctrl = tmp_ctrl(freq_ctr_idx,time_ctr_idx);
        tmp_ctrl = median(tmp_ctrl);
        tmp_exp = S_exp_downsampled{k};
        time_exp_idx = time_exp_dowmsampled{k} >= 15 & time_exp_dowmsampled{k} <=55;
        tmp_exp = tmp_exp(freq_exp_idx,time_exp_idx);
        tmp_exp = median(tmp_exp);
        S_ctr_to_test = [S_ctr_to_test;tmp_ctrl(:)];
        S_exp_to_test = [S_exp_to_test;tmp_exp(:)];
    end
    p_val_uncorrected_bands(i) = ranksum(S_ctr_to_test(:),S_exp_to_test(:));
    to_bar(1,i) = mean(S_ctr_to_test);
    to_bar(2,i) = mean(S_exp_to_test);
    err(1,i) = std(S_ctr_to_test)/sqrt(length(S_ctr_to_test));
    err(2,i) = std(S_exp_to_test)/sqrt(length(S_exp_to_test));
end
p_val_corrected_bands = bonf_holm(p_val_uncorrected_bands);

figure('units','normalized','outerposition',[0 0 1 1]);p = bar(to_bar');xticklabels(labels);
ngroups = 7;
nbars = 2;
% Calculating the width for each bar group
hold all
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = [];
for i = 1:nbars
    x(i,:) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x(i,:), to_bar(i,:), err(i,:), '.');
end
sigstar({[x(1,1),x(2,1)],[x(1,2),x(2,2)],[x(1,3),x(2,3)],[x(1,4),x(2,4)],[x(1,5),x(2,5)],[x(1,6),x(2,6)],[x(1,7),x(2,7)]},p_val_corrected_bands');
ax = gca;ax.FontSize = 30;
legend([p(1),p(2)],'Ctrl','Exp')

%% MUA statistics 
clear 
clc
close all

%% load results
path = ''; % add the path for the circular statistics toolbox
addpath(path);
[file,path] = uigetfile('*.mat','Select One or More results data', ...
    'MultiSelect', 'on'); % load the results that are obtained from the spiking_properties routine
Data = cell(1,length(file));
mousename = cell(1,length(file));
for i=1:length(file)
    Data{i} = load([path,file{i}]);
    mousename{i} = file{i}(1:end-4);
end

exp_idx = find(contains(mousename,'exp'));
ctrl_idx = find(contains(mousename,'ctr'));

conditions = {'baseline','transition','active'};

n_exp = length(exp_idx);
n_ctrl = length(ctrl_idx);
n_conditions = length(conditions);

%% channel lvl spike rate analysis
spk_rate_exp = [];
for cnd=1:n_conditions
    tmp_cnd = [];
    for id=1:n_exp
        temp = cell2mat(Data{exp_idx(id)}.spkrate.(conditions{cnd}))';
        tmp_cnd = cat(1,tmp_cnd,temp);
    end
    spk_rate_exp(:,cnd) = tmp_cnd;
end

spk_rate_ctrl = [];
for cnd=1:n_conditions
    tmp_cnd = [];
    for id=1:n_ctrl       
        temp = cell2mat(Data{ctrl_idx(id)}.spkrate.(conditions{cnd}))';
        tmp_cnd = cat(1,tmp_cnd,temp);       
    end
    spk_rate_ctrl(:,cnd) = tmp_cnd;
end

figure();
n_units_exp = 80; % total number of units in exp
n_units_ctrl = 80;% total number of units in ctrl
subplot(1,2,1);scatter(spk_rate_exp(:,1),spk_rate_exp(:,3));xlim([0,n_units_exp]);ylim([0,n_units_exp]);
hold all; plot(0:0.1:n_units_exp,0:0.1:n_units_exp);axis square;title('exp data');xlabel('Baseline spike rate');ylabel('Active spike rate');
subplot(1,2,2);scatter(spk_rate_ctrl(:,1),spk_rate_ctrl(:,3));xlim([0,n_units_ctrl]);ylim([0,n_units_ctrl]);
hold all; plot(0:0.1:n_units_ctrl,0:0.1:n_units_ctrl);axis square;title('ctrl data');xlabel('Baseline spike rate');ylabel('Active spike rate');

%% phase coupling analyses
count = 1;
phase_pool = [];
for id=1:n_exp
    temp = Data{exp_idx(id)}.phase;
    for cnd=1:n_conditions
        tmp_phase_pooled = temp.(conditions{cnd});
        PLV.exp.(conditions{cnd}){count} = cellfun(@circ_r,tmp_phase_pooled);
        mean_phase.exp.(conditions{cnd}){count} = cellfun(@circ_mean,tmp_phase_pooled);
        phase_pool.exp.(conditions{cnd}){count} = cell2mat(tmp_phase_pooled');
    end
    count = count + 1;
end

count = 1;
for id=1:n_ctrl
    temp = Data{ctrl_idx(id)}.phase;
    for cnd=1:n_conditions
        tmp_phase_pooled = temp.(conditions{cnd});
        PLV.ctrl.(conditions{cnd}){count} = cellfun(@circ_r,tmp_phase_pooled);
        mean_phase.ctrl.(conditions{cnd}){count} = cellfun(@circ_mean,tmp_phase_pooled);
        phase_pool.ctrl.(conditions{cnd}){count} = cell2mat(tmp_phase_pooled');       
    end
    count = count + 1;
end


PLV_gain_exp = (cell2mat(PLV.exp.active) - cell2mat(PLV.exp.baseline))./(cell2mat(PLV.exp.active) + cell2mat(PLV.exp.baseline));
PLV_gain_ctrl = (cell2mat(PLV.ctrl.active) - cell2mat(PLV.ctrl.baseline))./(cell2mat(PLV.ctrl.active) + cell2mat(PLV.ctrl.baseline));
spk_rate_gain_exp = (spk_rate_exp(:,3) - spk_rate_exp(:,1))./(spk_rate_exp(:,3) + spk_rate_exp(:,1));
spk_rate_gain_ctrl = (spk_rate_ctrl(:,3) - spk_rate_ctrl(:,1))./(spk_rate_ctrl(:,3) + spk_rate_ctrl(:,1));
figure();
subplot(1,2,1);scatter(spk_rate_gain_exp,PLV_gain_exp);ylim([-1,1]);xlim([-1,1]);axis square
title('exp data');xlabel('spike rate modulation index');ylabel('PLV modulation index');
subplot(1,2,2);scatter(spk_rate_gain_ctrl,PLV_gain_ctrl);ylim([-1,1]);xlim([-1,1]);axis square
title('ctrl data');xlabel('spike rate modulation index');ylabel('PLV modulation index');

disp('percentage of MUA units that decrease firing rate in exp data')
disp(length(find(spk_rate_gain_exp<0))./length(spk_rate_exp)*100)
disp('percentage of MUA units that decrease firing rate in ctrl data')
disp(length(find(spk_rate_gain_ctrl<0))./length(spk_rate_ctrl)*100)

figure()
edges = -pi/8:pi/4:2*pi+pi/8;
centers = movmean(edges,2);
centers = centers(2:end-1);
for cnd=1:n_conditions
    tmp_exp = phase_pool.exp.(conditions{cnd});
    to_hist = [];
    for i=1:length(tmp_exp)
        to_hist = [to_hist;tmp_exp{i}];
    end
    to_hist = wrapTo2Pi(to_hist);
    [counts,edges] = histcounts(to_hist,edges);
    counts(1) = counts(1)+counts(end);
    counts(end) = [];
    prob = counts./sum(counts);
    subplot(2,3,cnd);bar(rad2deg(centers),prob);xticks(rad2deg(centers));
    title(conditions{cnd});ylim([0,0.3]);ylabel('probability of spike occurance');
    xlabel('Angle of LFP (deg)'); axis square
    
    to_hist = [];
    tmp_ctrl = phase_pool.ctrl.(conditions{cnd});
    for i=1:length(tmp_ctrl)
        to_hist = [to_hist;tmp_ctrl{i}];
    end
    to_hist = wrapTo2Pi(to_hist);
    [counts,edges] = histcounts(to_hist,edges);
    counts(1) = counts(1)+counts(end);
    counts(end) = [];
    prob = counts./sum(counts);
    subplot(2,3,cnd+3);bar(rad2deg(centers),prob);xticks(rad2deg(centers));
    ylabel('probability of spike occurance');xlabel('Angle of LFP (deg)');ylim([0,0.3]);axis square
end


to_bar = [mean(PLV_gain_exp),mean(PLV_gain_ctrl)];
err_bar = [std(PLV_gain_exp)/length(PLV_gain_exp),std(PLV_gain_ctrl)/length(PLV_gain_ctrl)];
p_val_PLV = ranksum(PLV_gain_exp,PLV_gain_ctrl);
figure();bar(to_bar);hold all;errorbar(to_bar,err_bar,'.');sigstar([1,2],p_val_PLV)
xticklabels({'exp','ctrl'});ylabel('PLV modulation index');

m_active_exp = cell2mat(mean_phase.exp.active)';
r_active_exp = cell2mat(PLV.exp.active)';
disp('percentage of MUA units with significant coupling in exp data')
disp(sum(r_active_exp>0.1)./length(spk_rate_exp)*100)
figure();polarhistogram(m_active_exp(r_active_exp>0.1),10);title('preferred phase of MUA units with non-uniform distribution exp data')


m_active_ctrl = cell2mat(mean_phase.ctrl.active)';
r_active_ctrl = cell2mat(PLV.ctrl.active)';
disp('percentage of MUA units with significant coupling in ctrl data')
disp(sum(r_active_ctrl>0.1)./length(spk_rate_ctrl)*100)
figure();polarhistogram(m_active_ctrl(r_active_ctrl>0.1),10);title('preferred phase of MUA units with non-uniform distribution ctrl data')

figure()
t = linspace(0,2*pi,10000);
plot(rad2deg(t),cos(t),'LineWidth',2)
hold all
tmp_to_plot = m_active_exp(r_active_exp>0.1);
plot(rad2deg(wrapTo2Pi(tmp_to_plot)),cos(tmp_to_plot)+0.1,'o','LineWidth',2)
tmp_to_plot = m_active_ctrl(r_active_ctrl>0.1);
plot(rad2deg(wrapTo2Pi(tmp_to_plot)),cos(tmp_to_plot)-0.1,'o','LineWidth',2)
xticks(rad2deg(centers));axis square;xlabel('Angle (deg)');yticks([])