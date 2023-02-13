function allExptData = loadLatestComboSensorData()
    % CD to the combo sensor experiment directory before this
    file_data = dir(fullfile('**', 'all_combo_expt_data*.mat'));
    

    if isempty(file_data)
        disp('No combo sensor data found.');
    elseif length(file_data)==1
        load(fullfile(file_data(1).folder, file_data(1).name), 'allExptData');
    elseif length(file_data)>1
        tempnames = {file_data.name};
        tempdates = NaT(1, length(tempnames));
        for i = 1:length(tempnames)
            tempdat = erase(tempnames{i}, 'all_combo_expt_data');
            tempdat = erase(tempdat, '.mat');
            if ~isempty(tempdat)
                tempdat = tempdat(2:end); % get rid of leading _
                tempdates(i) = datetime(tempdat);
            end
        end
        % get the latest date
        [~, ind] = max(tempdates);
        % load that one
        load(fullfile(file_data(ind).folder, file_data(ind).name), 'allExptData');
    end
            



end