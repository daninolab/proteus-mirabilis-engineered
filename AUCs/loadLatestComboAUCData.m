function comboGeneAUCData = loadLatestComboAUCData()
    % CD to the timelapse directory
    file_data = dir(fullfile('**', 'combo_AUCs_*.mat'));

    if isempty(file_data)
        disp('No combo gene experiment data found.');
    elseif length(file_data)==1
        load(fullfile(file_data(1).folder, file_data(1).name), 'comboGeneAUCData');
    elseif length(file_data)>1
        tempnames = {file_data.name};
        tempdates = NaT(1, length(tempnames));
        for i = 1:length(tempnames)
            tempdat = erase(tempnames{i}, 'combo_AUCs_');
            tempdat = erase(tempdat, '.mat');
            if ~isempty(tempdat)
                if strcmpi(tempdat(1), '_')
                    tempdat = tempdat(2:end); % get rid of leading _
                end
                tempdates(i) = datetime(tempdat);
            end
        end
        % get the latest date
        [~, ind] = max(tempdates);
        % load that one
        load(fullfile(file_data(ind).folder, file_data(ind).name), 'comboGeneAUCData');
    end
            



end