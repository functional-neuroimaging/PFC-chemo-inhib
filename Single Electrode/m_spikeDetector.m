function [spkTimes, spkIndex, spkShape, threshold] = m_spikeDetector(MUA, fsMUA, method, th)
% M_SPIKEDETECTOR Detect and extract spike waveform from MUA signal.
% Detection method is based on threshold which can be given from the
% experimentalist or predefined by the script itself.
% DEFAULT TH: 'median'. Another option: 'std'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline:
% - Define th
% - Find spikes over th
% - Find the highest peak of the spike
% - Export the spikes shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   MUA : Multi unit activity
%   fsMUA : (optional) MUA sampling frequency. DEFAULT: 4000 Hz 
%   method : (optional) how to determine the threshold. 
%               option: 'median', '5std'. DEFAULT: 'median'
%   th : (optional) threshold to utilize for the detection 
%   
% OUTPUT:
%   spkTimes : spike position in second
%   spkIndex : spike position in index
%   spkShape : spike shape
%   threshold : used threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changelog:
%   v1.1 24/06/2019 FR: Bug about external th corrected
%   v1.0 15/05/2019 FR: function created 

% checking input args
if nargin < 1
    error('Missing input arguments');
elseif nargin < 2
    fsMUA = 4000;
    ext_th = false;
    method = 'median';
elseif nargin < 3
    ext_th = false;
    method = 'median';
elseif nargin < 4 
    ext_th = false;
else
    ext_th = true; 
end

%init
threshold = zeros(size(MUA,2),1);

for ch = 1:size(MUA,2)
    % threshold initialization
    if ext_th
        threshold = th;
    else
        switch method
            case 'median'
                coeff = 4;
                threshold(ch) = coeff*median(abs(MUA(:,ch))/0.6745);
            case '5std'
                coeff = 5;
                threshold(ch) = coeff*std(MUA(:,ch));
            otherwise
                error('Wrong spike detection method!')
        end
    end

    % find point over the th
    spike_index = find (abs(MUA(:,ch)) > threshold(ch));
    spike_width = 0.001;
    l = floor(spike_width * fsMUA);

    % finding only the largest peak
    if l > 0 
        oldlength = 0;
        while length(spike_index) ~= oldlength
            oldlength = length(spike_index);
            k = find(diff(spike_index)<l);
            j = find(abs(MUA(spike_index(k),ch)) <= abs(MUA(spike_index(k+1),ch)));
            spike_index(k(j) )= [];
            k = find(diff(spike_index)<l);
            j = find(abs(MUA(spike_index(k+1),ch)) < abs(MUA(spike_index(k),ch)));
            spike_index(k(j)+1) = [];
        end
    end

    spike_index(spike_index < l+1) = [];
    spike_index(spike_index > size(MUA,1) - l-1) = [];
    spike_times = (spike_index - 1)/fsMUA;

    %% spikes shape 
    spikes_shape = zeros(l+1,length(spike_index));
    for j = 1:length(spike_index)
        spikes_shape(:,j) = MUA(spike_index(j)-l/2:spike_index(j)+l/2,ch);
    end

    spkTimes{ch} = spike_times;
    spkIndex{ch} = spike_index;
    spkShape{ch} = spikes_shape;
end
end

