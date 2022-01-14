function [data, mouseName, path] = io_importProcessedData
%IO_IMPORTPROCESSEDDATA function to import data selected by a UI. This
%function will import all the data inside a single structure with a generic
%name. The specific name is keep insiede mouseName with the same index.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%   data : struct containing a struct for each mouse imported. Index names are generic. 
%   mouseName : reference to real file name
%   path : path selected in the UI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changelog
%   v1.0 20/05/2019 FR: Function created

[file,path] = uigetfile('*.mat','Select One or More preprocessed ePHYS data', ...
    'MultiSelect', 'on');

if ~iscell(file)
    file = {file};
end
    
mouseName = cell(1,size(file,2));

for j = 1:size(file,2)
     tempData = load([path,file{j}]);
     mouseName(j) = fieldnames(tempData);
     data.(['d',num2str(j)]) = tempData.(mouseName{j});
     fprintf('%s imported!\n',file{j})
end
end


