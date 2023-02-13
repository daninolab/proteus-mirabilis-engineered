function updateTimelapseData(timelapsedir)
    % This function should check the given folder for a data structure with
    % timelapse analysis, then look for new timelapse analyses and add them
    % in to the structure. Use only the basic fields (can add new fields of
    % analysis to structure later/in other functions)

    % Check the folder for a file called "all_timelapse_data"-date-".mat"
    oldFolder = cd(timelapsedir);
    % Check the folder for a file called "all_timelapse_data"-date-".mat"
    file_data = dir(fullfile('**', 'all_timelapse_data*.mat'));
    
    if isempty(file_data)
        disp('Compiled timelapse data file does not exist');
        % return some value to main script?
        return;
    
    else
        % LOAD THE LATEST ALL TIMELAPSE DATA STRUCTURE
        disp('Found folder.');
%         disp(file_data(1).name);
        allExptData = loadLatestTimelapseData();
%         load(fullfile(file_data(1).folder, file_data(1).name), 'allExptData'); %allExptData

        % Look for all "scanlapse_2" files in the subdirectories
        file_data = dir(fullfile('**', 'scanlapse_msmts_2_*.mat'));
        file_names = {file_data.name}; %cell array of file names
        new_expts = false;
        
        % Load the dates of previously analyzed timelapses
        prev_dates = {allExptData.date};
        for i = 1:length(prev_dates)
            prev_dates{i} = formatDateStr(prev_dates{i});
        end
        
        % Check if any are not yet added to main structure, then load those and
        % add their data into main structure
        for i = 1:length(file_names)
            current = file_names{i};
            curr_date = erase(current, 'scanlapse_msmts_2_');
            curr_date = erase(curr_date, '.mat'); 
            curr_date = formatDateStr(curr_date);

            % check if any dates match it
            if ~any(cellfun(@(x) isequal(x, curr_date), prev_dates))
                % This experiment's analysis is not yet in data structure
                new_expts = true;
                disp(strcat("Found a new timelapse to add ", string(curr_date)));
                mainFolder = cd(file_data(i).folder);
                % ADD THE NEW ONE
                temp = getTimelapseData(current);
                % Add to the end of current struct
                k = length(allExptData)+1;
                allExptData(k).files = file_names{i};
                allExptData = addNewData(allExptData, temp, file_data(i).folder, k);
                cd(mainFolder);                
            end
                

        end
        if ~new_expts
            disp('No new expts found.');
        end
    end

    disp('Done. Saving');
    % save updated file data
    filename = strcat('all_timelapse_data_', date, '.mat');
%     filename = 'all_timelapse_data.mat';
    save(filename, '-v7.3');
    cd(oldFolder);
end



function allExptData = addNewData(allExptData, temp, foldname, k)
   disp('Adding basic data...');
   % Add a new timelapse analysis to current data structure
   allExptData(k).expt =temp{1};
   allExptData(k).date = temp{2};
   allExptData(k).exptfolder = foldname;
   
   disp('Adding file & plate data...');
   allExptData(k).num_files = temp{3};
   allExptData(k).plates = temp{4};
   allExptData(k).imfiles = temp{5};

   % Get the times of the data
   disp('Adding basic analysis data...');
   allExptData(k).times = makeTimeVec(foldname);
   allExptData(k).dists = temp{6};
   allExptData(k).avgd = temp{7};
   allExptData(k).cvs = temp{8};
   allExptData(k).stdevs = temp{9};
        
end

function times = makeTimeVec(foldname)
    % Cd to the folder of expt. Check if there is a timestamps.txt file. if
    % not, try getting all the scan times. If those don't look right, make
    % artificial vector.
    % Go to expt folder:
    oldFolder = cd(foldname);
    % Check for timestamps.txt:
    temp = dir(fullfile('**', 'timestamps.txt'));
    if ~isempty(temp)
        % Go to folder where timestamps.txt is
        cd(temp(1).folder);
        times = getTimesFromTxt();
        cd(foldname);
    else
        % First try getting path to scans & times of creation of those
        disp(foldname);
        tifpath = input('Input the path to the tif scans for expt: ', 's');
        times = getTimes(tifpath);
        if mean(diff(times))<0.01
            % this is saying images were taken or saved less than a minute
            % apart--probably original creation date not saved
            % generate a time vector artificially
            times = makeTimes(tifpath);
        end
    end
    
      
    cd(oldFolder);
end