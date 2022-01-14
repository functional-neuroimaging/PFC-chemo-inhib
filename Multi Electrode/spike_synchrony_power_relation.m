% routine to do the simulation to show coupling to specific frequency would
% increase the power, supplementary figure 7

clear 
clc
close all

%% LFP generation
fs = 1000;
t = 0:1/fs:10;
f = 1;
cos_L = cos(2*pi*f*t);
LFP_L = cos_L + 0.1*randn(size(cos_L));
angle_L = angle(hilbert(LFP_L));
figure();
subplot(3,2,1)
plot(t,cos_L,'LineWidth',1,'Color','k');
xlabel('time (seconds)')
hold all
%% spike generation - coupled case
rate = 10;
n_spikes = rate*t(end);
spike_count = zeros(1,length(angle_L));
coupling_strenght = 20;
preferred_phase = -pi;

alpha = circ_vmrnd(preferred_phase, coupling_strenght, n_spikes);
idx = knnsearch(angle_L',alpha,'K',10);
perm = zeros(size(idx));
for i=1:length(idx)
    idx_perm = randperm(size(idx,2),1);
    perm(i,idx_perm) = 1;
end
idx_spike = idx(logical(perm));
spike_count(idx_spike) = 1;
smoothing_window = gausswin(500);
spike_train = filter(smoothing_window,1,spike_count);
subplot(3,2,5);
spk = t(idx_spike);
plot([spk; spk], repmat(ylim',1,size(spk,2)), '-b')
ylim([0,2]);
subplot(1,2,2)
pspectrum(spike_train - mean(spike_train),fs,'FrequencyLimits',[0,5]);
%% spike generation - uncoupled case
rate = 10;
n_spikes = rate*t(end);
%spike_count = zeros(1,length(angle_L));
coupling_strenght = 0.5;
preferred_phase = 0;
alpha = circ_vmrnd(preferred_phase, coupling_strenght, n_spikes);
idx = knnsearch(angle_L',alpha,'K',10);
perm = zeros(size(idx));
for i=1:length(idx)
    idx_perm = randperm(size(idx,2),1);
    perm(i,idx_perm) = 1;
end
idx_spike = idx(logical(perm));
spike_count(idx_spike) = 1;
idx_spike_all = find(spike_count);
smoothing_window = gausswin(500);
spike_train = filter(smoothing_window,1,spike_count);
subplot(3,2,3);
spk = t(idx_spike_all);
plot([spk; spk], repmat(ylim',1,size(spk,2)), '-r')
ylim([0,2]);
subplot(1,2,2)
hold all
pspectrum(spike_train - mean(spike_train),fs,'FrequencyLimits',[0,5]);


