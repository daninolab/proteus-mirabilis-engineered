function temp = getSingleExptData(curr_file)
    % Given the filename of a '...radial_analysis.mat' file, load the
    % relevant data into temp
    temp = {};

    % First, check if the file has polar_analysis (new)
    listOfVariables = who('-file', curr_file);

    if ismember('polar_analysis', listOfVariables)
    
        % First, load polar_analysis, which is in both 'analysis_polar.mat' and
        % 'radial_analysis.mat' files
        disp('Loading expt data...');        
        load(curr_file, 'polar_analysis');
        disp('loaded');

        % Polar analysis struct will have mostly same fields but a few will
        % differ if old or new. Old one has 'useim' field but empty everywhere
        % except in that field for the images not to be used, so instead we can
        % look for empty indices in the filename
        % HOWEVER some analysis was done in the v v old format with different
        % structs 's', 'a', 'analysis'. Check for those just as in the combo
        % sensor expts
        picinds = find(~cellfun('isempty', {polar_analysis.filename}));
        temp{1} = {polar_analysis(picinds).filename};
        temp{3} = [polar_analysis(picinds).colcenter];

        % The colony radius in pixels is only given in the new analysis ones.
        % Check if it is a field, otherwise keep blank
        if isfield(polar_analysis, 'colrad')
            temp{4} = [polar_analysis(picinds).colrad];
        else
            disp('Colony radius not analyzed in this experiment, leaving empty.');
        end

        temp{5} = {polar_analysis(picinds).petrimask};
        temp{6} = {polar_analysis(picinds).dist_vec};
        temp{7} = {polar_analysis(picinds).avgd};
        temp{8} = {polar_analysis(picinds).cvs};
        temp{9} = {polar_analysis(picinds).stdevs};
    elseif ismember('s', listOfVariables) 
        % This is old old style analysis
        load(curr_file, 's');
        load(curr_file, 'analysis');
        picinds = find([s.usepic]==1);
        
        % Start filling in the temp cell
        temp{1} = {s(picinds).img_name};
        temp{3} = [s(picinds).col_center];
        % No petri mask saved--leave temp{5} blank        
        temp{6} = {s(picinds).dists};
        temp{7} = {s(picinds).avgd};
        temp{8} = {analysis(picinds).cvs};
        temp{9} = {analysis(picinds).stdevs};
        
    end
    
    % For temp 2, we want the polar image names, these can be created from
    % the original scan image names in temp{1}

    for k = 1:length(temp{1})
        tempfile = temp{1}{k};
        tempname = erase(tempfile, '.tif');
        newname = strcat(tempname, '_polarim.tif');
        temp{2}{end+1} = newname;
    end
end