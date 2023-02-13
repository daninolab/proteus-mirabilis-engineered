function updateSingleGeneData(singleExptDir, savewhendone)
    % This function should check the given folder for a data structure with
    % single gene expt analysis, then look for new analyses and add them
    % in to the structure. Use only the basic fields (can add new fields of
    % analysis to structure later/in other functions)

    oldFolder = cd(singleExptDir);
    % Check the folder for a file called "all_timelapse_data"-date-".mat"
    file_data = dir(fullfile('**', 'all_single_gene_expt_data*.mat'));
    if isempty(file_data)
            disp('Compiled timelapse data file does not exist');
            % return some value to main script?
            return;
    else
        count = 0;
        allExptData = loadLatestSingleGeneData();
        
        % Now that we have loaded it, let's check for analyses that have
        % yet to be added
        % Look for all "scanlapse_2" files in the subdirectories
        file_data = dir(fullfile('**', '*radial_analysis.mat'));
        file_names = {file_data.name}; %cell array of file names
        file_data = dir(fullfile('**', '*analysis_polar.mat'));
        file_names = [file_names, {file_data.name}];
        file_data = dir(fullfile('**', fullfile('**', '*analysis_polar.mat')));
        file_names = [file_names, {file_data.name}];
        file_names = unique(file_names, 'stable');
        new_expts = false;

        % Check the exptfiles contained in allExptData

        done_files = {allExptData.exptfile};
        for i = 1:length(file_names)
            if any(contains(done_files, file_names{i}))
                disp('found');
                disp(file_names{i});
            else
                disp('not found');
                disp(file_names{i});
                % now decide whether to add this or not
            end
        end
        
        
%         prev_dates = {allExptData.date};
%         for i = 1:length(prev_dates)
%             prev_dates{i} = formatDateStr(prev_dates{i});
%         end
        
         % Check if any are not yet added to main structure, then load those and
        % add their data into main structure
        for i = 1:length(file_names)
            current = file_names{i};
            if any(contains(done_files, current))
                disp('found');
                disp(current);
            else            
                % This experiment's analysis is not yet in data structure
                new_expts = true;
                disp('Found an expt not added yet: ');
                disp(current);
                resp = input('Add this expt? Type y/n: ', 's');
                
                if strcmpi(resp, 'y') | strcmpi(resp, 'yes')
                    count = count+1;
                    % Get the folder it is in
                    file_data = dir(fullfile('**', file_names{i}));
                    mainFolder = cd(file_data(1).folder);
                    % ADD THE NEW ONE
                    temp = getSingleExptData(current);
                    % Add to the end of current struct
                    k = length(allExptData)+1;
                    allExptData(k).exptfile = file_names{i};
                    allExptData = addNewData(allExptData, temp, file_data(1).folder, k);
                    cd(mainFolder);  
                end
            end

            
        end
        if ~new_expts
            disp('No new expts found.');
        end
        fprintf('Added %g experiments.\n', count);
    end
    
    if savewhendone
        disp('Done. Saving...');
        % save updated file data
        filename = strcat('all_single_gene_expt_data_', date, '.mat');
        save(filename, 'allExptData', '-v7.3');
    else
        disp('Done. Not saving.');
    end
end



function allExptData = addNewData(allExptData, temp, foldname, k)
    disp('Adding basic data...');
    % Add a new single gene analysis to current data structure

    allExptData(k).exptfolder = foldname;
    disp(foldname);
    allExptData(k).date = formatDateStr(input('Type in expt date: ', 's'));
    % Get the scanfolder
    tempfolders = split(foldname, '/');
    allExptData(k).scanfolder = fullfile(tempfolders{1:end-1});
    
    disp('Adding file & plate data...');

    allExptData(k).imfiles = temp{1}; % store tif file names
    allExptData(k).polfiles = temp{2}; % store names of polar files
    allExptData(k).colcenters = temp{3};
    allExptData(k).colrads_pix = temp{4};
    allExptData(k).petrimasks = temp{5};
    allExptData(k).dists = temp{6};
    allExptData(k).avgd = temp{7};
    allExptData(k).cvs = temp{8};
    allExptData(k).stdevs = temp{9};
   
        
end