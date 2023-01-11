function [sampledData, sampledVicon, sampledTime,proj2_Data] = init(dataNum)
%RESAMPLEDATA Resample the given Vicon Data
data1 = 'studentdata1.mat';
data4 = 'studentdata4.mat';

switch dataNum
    case 1
        allData = load(data1);
        proj2_Data = load('proj2_dataset1.mat');
    case 4
        allData = load(data4);
        proj2_Data = load('proj2_dataset4.mat');
end
        
dataTimes = vertcat(allData.data(:).t); % Extract all values of t field from data
matchingIndices = ismember(allData.time, dataTimes); % Obtain incdices of matches
sampledTime = allData.time(matchingIndices); % Keep only 1 indices
sampledData = allData.data;
sampledVicon = allData.vicon(:, matchingIndices); % Keep only 1 indices

end