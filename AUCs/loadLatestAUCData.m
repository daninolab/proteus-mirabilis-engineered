function singleGeneAUCData = loadLatestAUCData()
    % CD to the timelapse directory
    file_data = dir(fullfile('**', 'single_gene_AUCs_*.mat'));
    

    if isempty(file_data)
        disp('No single gene experiment data found.');
    elseif length(file_data)==1
        load(fullfile(file_data(1).folder, file_data(1).name), 'singleGeneAUCData');
    elseif length(file_data)>1
        tempnames = {file_data.name};
        tempdates = NaT(1, length(tempnames));
        for i = 1:length(tempnames)
            tempdat = erase(tempnames{i}, 'single_gene_AUCs_');
            tempdat = erase(tempdat, '.mat');
            if ~isempty(tempdat)
                tempdat = tempdat(2:end); % get rid of leading _
                tempdates(i) = datetime(tempdat);
            end
        end
        % get the latest date
        [~, ind] = max(tempdates);
        % load that one
        load(fullfile(file_data(ind).folder, file_data(ind).name), 'singleGeneAUCData');
    end
            



end