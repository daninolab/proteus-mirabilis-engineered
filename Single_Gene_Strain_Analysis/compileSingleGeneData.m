function compileSingleGeneData(singlegene_exptdir)
    % This function should check the given folder for a data structure with
    % single gene analysis, if it does not exist, it should create one. In
    % the data structure, it should fill with available single gene expt
    % analysis from different '...radial_analysis.mat' OR
    % '...analysis_polar.mat' files. Since those are probably in a
    % different format, be ready to format....
    
    % Want to make a file 'all_single_expt_data_{date}.mat' at the end
    
    oldFolder = cd(singlegene_exptdir);
    % Check the folder for a file called "all_timelapse_data"-date-".mat"
    file_data = dir(fullfile('**', 'all_single_expt_data_{date}.mat'));
    
    if ~isempty(file_data)
        disp('Compiled combo data file already exists');
        % return some value to main script?
        return;
    else
        % If it does not exist, create new structure
        disp('Creating new data structure...');
        allExptData = struct;

        % Look for all analysis files in the subdirectories
        file_data1 = dir(fullfile('**', '*radial_analysis.mat'));
        
        % Look for older style analysis files
        file_data2 = dir(fullfile('**', '*analysis_polar.mat'));
        
        % Combine
        file_data = [file_data1; file_data2];
        file_names = {file_data.name}; %cell array of file names
        fprintf('Found %d single gene expt files \n', length(file_names));
        
        % Load each and use the data to fill in the main structure
        for i = 1:length(file_names)
           disp(strcat("Adding data from: ", file_names{i}));
           mainFolder = cd(file_data(i).folder);
           
           %Fill in basic data for that expt
           allExptData(i).exptfile = file_names{i};
           allExptData(i).exptfolder = file_data(i).folder;
           % Look for date
           tempstrs = strsplit(file_names{i}, '_');
           if length(tempstrs)==4 && strcmpi(tempstrs{3}, 'radial') && strcmpi(tempstrs{4}, 'analysis.mat')
               tempdate = formatDateStr(tempstrs{1});
           else
               tempdate = formatDateStr(input('Type in expt date: ', 's'));
           end
           allExptData(i).date = tempdate;
           
           % Load the major analysis
           current = file_names{i};
           temp = getSingleExptData(current);
           
           % Fill in allExptData using Temp Struct

           allExptData(i).imfiles = temp{1}; % store tif file names
           allExptData(i).polfiles = temp{2}; % store names of polar files
           allExptData(i).colcenters = temp{3};
           allExptData(i).colrads_pix = temp{4};
           allExptData(i).petrimasks = temp{5};
           allExptData(i).dists = temp{6};
           allExptData(i).avgd = temp{7};
           allExptData(i).cvs = temp{8};
           allExptData(i).stdevs = temp{9};

           
           cd(mainFolder);

        end % End iterating over all the exptfiles

        disp('Completed compiling timelapse data. Saving to current folder...');
        % SAVE THE STRUCT TO A FILE
        % Save the measurements in a .mat file in corresponding folder, use a
        % version compatible with large variables. Note: future updates should
        % be saved with the date in the name.
        filename = 'all_single_gene_expt_data.mat';
        save(filename, 'allExptData', '-v7.3');
        cd(oldFolder);
        disp('Done.');

    end % End 'if isempty filedata'
    
end



