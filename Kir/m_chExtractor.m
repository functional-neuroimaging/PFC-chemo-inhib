function [splittedStruct] = m_chExtractor(preprocessedStruct, chToKeep)
%M_CHEXTRACTOR Extract from a preprocessed struct only specified channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   preprocessedStruct : preprocessed struct of data
%   chToKeep : channel to keep
%   
% OUTPUT:
%   splittedStruct : new struct with only the channel given in chToKeep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changelog:
%   v1.0 24/06/2019 FR: Function created

%raw
preprocessedStruct.raw.data = preprocessedStruct.raw.data(:, chToKeep);
%LFP
preprocessedStruct.LFP.data = preprocessedStruct.LFP.data(:, chToKeep);
preprocessedStruct.LFP.dataNo50Hz = preprocessedStruct.LFP.dataNo50Hz(:, chToKeep);
%MUA
preprocessedStruct.MUA.data = preprocessedStruct.MUA.data(:, chToKeep);
preprocessedStruct.MUA.analg_MUA = preprocessedStruct.MUA.analg_MUA(:, chToKeep);
preprocessedStruct.MUA.unfilt_MUA = preprocessedStruct.MUA.unfilt_MUA(:, chToKeep);
%spikes
preprocessedStruct.spikes = preprocessedStruct.spikes(:, chToKeep);
%output
splittedStruct = preprocessedStruct;
end

