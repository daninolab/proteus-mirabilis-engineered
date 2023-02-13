%% For the Revision, Plot & Format copA/cadA Plate Reader Data

%% Load the OD data
opts = spreadsheetImportOptions("NumVariables", 97);

% Specify sheet and range
opts.Sheet = "OD Data Analysis For Matlab";
opts.DataRange = "BE8:EW59";

% Specify column names and types
opts.VariableNames = ["Timeh1", "LB", "LB1", "LB01UMCuSO4", "LB01UMCuSO1", "LB1UMCuSO4", "LB1UMCuSO1", "LB10UMCuSO4", "LB10UMCuSO1", "LB100UMCuSO4", "LB100UMCuSO1", "LB1MMCuSO4", "LB1MMCuSO1", "PM", "PM1", "PM01UMCuSO4", "PM01UMCuSO1", "PM1UMCuSO4", "PM1UMCuSO1", "PM10UMCuSO4", "PM10UMCuSO1", "PM100UMCuSO4", "PM100UMCuSO1", "PM1MMCuSO4", "PM1MMCuSO1", "PMCopAGFP", "PMCopAGFP1", "PMCopAGFP01UMCuSO4", "PMCopAGFP01UMCuSO1", "PMCopAGFP1UMCuSO4", "PMCopAGFP1UMCuSO1", "PMCopAGFP10UMCuSO4", "PMCopAGFP10UMCuSO1", "PMCopAGFP100UMCuSO4", "PMCopAGFP100UMCuSO1", "PMCopAGFP1MMCuSO4", "PMCopAGFP1MMCuSO1", "PM2", "PM3", "PM01UMZnCl2", "PM01UMZnCl1", "PM1UMZnCl2", "PM1UMZnCl1", "PM10UMZnCl2", "PM10UMZnCl1", "PM100UMZnCl2", "PM100UMZnCl1", "PM1MMZnCl2", "PM1MMZnCl1", "PMCadAGFP", "PMCadAGFP1", "PMCadAGFP01UMZnCl2", "PMCadAGFP01UMZnCl1", "PMCadAGFP1UMZnCl2", "PMCadAGFP1UMZnCl1", "PMCadAGFP10UMZnCl2", "PMCadAGFP10UMZnCl1", "PMCadAGFP100UMZnCl2", "PMCadAGFP100UMZnCl1", "PMCadAGFP1MMZnCl2", "PMCadAGFP1MMZnCl1", "PMCadAGFP2", "PMCadAGFP3", "PMCadAGFP01UMZnCl3", "PMCadAGFP01UMZnCl4", "PMCadAGFP1UMZnCl3", "PMCadAGFP1UMZnCl4", "PMCadAGFP10UMZnCl3", "PMCadAGFP10UMZnCl4", "PMCadAGFP100UMZnCl3", "PMCadAGFP100UMZnCl4", "PMCadAGFP1MMZnCl3", "PMCadAGFP1MMZnCl4", "PM4", "PM5", "PM01UMHgCl2", "PM01UMHgCl1", "PM1UMHgCl2", "PM1UMHgCl1", "PM10UMHgCl2", "PM10UMHgCl1", "PM100UMHgCl2", "PM100UMHgCl1", "PMCadAGFP4", "PMCadAGFP5", "PMCadAGFP01UMHgCl2", "PMCadAGFP01UMHgCl1", "PMCadAGFP1UMHgCl2", "PMCadAGFP1UMHgCl1", "PMCadAGFP10UMHgCl2", "PMCadAGFP10UMHgCl1", "PMCadAGFP100UMHgCl2", "PMCadAGFP100UMHgCl1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
copAcadA_ODtable = readtable("/Users/anjali/Dropbox/_Shared_Patterning/_final_swarm_code_ncb/Datasets_for_Upload/10-02-22_copA-cadA_upload.xlsx", opts, "UseExcel", false);

clear opts;

%% Load the GFP data

opts = spreadsheetImportOptions("NumVariables", 97);

% Specify sheet and range
opts.Sheet = "OD Data Analysis For Matlab";
opts.DataRange = "BE304:EW355";

% Specify column names and types
opts.VariableNames = ["Timeh1", "LB", "LB1", "LB01UMCuSO4", "LB01UMCuSO1", "LB1UMCuSO4", "LB1UMCuSO1", "LB10UMCuSO4", "LB10UMCuSO1", "LB100UMCuSO4", "LB100UMCuSO1", "LB1MMCuSO4", "LB1MMCuSO1", "PM", "PM1", "PM01UMCuSO4", "PM01UMCuSO1", "PM1UMCuSO4", "PM1UMCuSO1", "PM10UMCuSO4", "PM10UMCuSO1", "PM100UMCuSO4", "PM100UMCuSO1", "PM1MMCuSO4", "PM1MMCuSO1", "PMCopAGFP", "PMCopAGFP1", "PMCopAGFP01UMCuSO4", "PMCopAGFP01UMCuSO1", "PMCopAGFP1UMCuSO4", "PMCopAGFP1UMCuSO1", "PMCopAGFP10UMCuSO4", "PMCopAGFP10UMCuSO1", "PMCopAGFP100UMCuSO4", "PMCopAGFP100UMCuSO1", "PMCopAGFP1MMCuSO4", "PMCopAGFP1MMCuSO1", "PM2", "PM3", "PM01UMZnCl2", "PM01UMZnCl1", "PM1UMZnCl2", "PM1UMZnCl1", "PM10UMZnCl2", "PM10UMZnCl1", "PM100UMZnCl2", "PM100UMZnCl1", "PM1MMZnCl2", "PM1MMZnCl1", "PMCadAGFP", "PMCadAGFP1", "PMCadAGFP01UMZnCl2", "PMCadAGFP01UMZnCl1", "PMCadAGFP1UMZnCl2", "PMCadAGFP1UMZnCl1", "PMCadAGFP10UMZnCl2", "PMCadAGFP10UMZnCl1", "PMCadAGFP100UMZnCl2", "PMCadAGFP100UMZnCl1", "PMCadAGFP1MMZnCl2", "PMCadAGFP1MMZnCl1", "PMCadAGFP2", "PMCadAGFP3", "PMCadAGFP01UMZnCl3", "PMCadAGFP01UMZnCl4", "PMCadAGFP1UMZnCl3", "PMCadAGFP1UMZnCl4", "PMCadAGFP10UMZnCl3", "PMCadAGFP10UMZnCl4", "PMCadAGFP100UMZnCl3", "PMCadAGFP100UMZnCl4", "PMCadAGFP1MMZnCl3", "PMCadAGFP1MMZnCl4", "PM4", "PM5", "PM01UMHgCl2", "PM01UMHgCl1", "PM1UMHgCl2", "PM1UMHgCl1", "PM10UMHgCl2", "PM10UMHgCl1", "PM100UMHgCl2", "PM100UMHgCl1", "PMCadAGFP4", "PMCadAGFP5", "PMCadAGFP01UMHgCl2", "PMCadAGFP01UMHgCl1", "PMCadAGFP1UMHgCl2", "PMCadAGFP1UMHgCl1", "PMCadAGFP10UMHgCl2", "PMCadAGFP10UMHgCl1", "PMCadAGFP100UMHgCl2", "PMCadAGFP100UMHgCl1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
copAcadA_GFPtable = readtable("/Users/anjali/Dropbox/_Shared_Patterning/_final_swarm_code_ncb/Datasets_for_Upload/10-02-22_copA-cadA_upload.xlsx", opts, "UseExcel", false);

clear opts;

%% Get the variable names and concentrations

all_strain_info_order = {};
% Let's make row 1 of this be the strain
% Row 2 will be the promoter
% row 3 will be the metal
all_vars = ["Timeh1", "LB", "LB1", "LB01UMCuSO4", "LB01UMCuSO1", "LB1UMCuSO4", "LB1UMCuSO1", "LB10UMCuSO4", "LB10UMCuSO1", "LB100UMCuSO4", "LB100UMCuSO1", "LB1MMCuSO4", "LB1MMCuSO1", "PM", "PM1", "PM01UMCuSO4", "PM01UMCuSO1", "PM1UMCuSO4", "PM1UMCuSO1", "PM10UMCuSO4", "PM10UMCuSO1", "PM100UMCuSO4", "PM100UMCuSO1", "PM1MMCuSO4", "PM1MMCuSO1", "PMCopAGFP", "PMCopAGFP1", "PMCopAGFP01UMCuSO4", "PMCopAGFP01UMCuSO1", "PMCopAGFP1UMCuSO4", "PMCopAGFP1UMCuSO1", "PMCopAGFP10UMCuSO4", "PMCopAGFP10UMCuSO1", "PMCopAGFP100UMCuSO4", "PMCopAGFP100UMCuSO1", "PMCopAGFP1MMCuSO4", "PMCopAGFP1MMCuSO1", "PM2", "PM3", "PM01UMZnCl2", "PM01UMZnCl1", "PM1UMZnCl2", "PM1UMZnCl1", "PM10UMZnCl2", "PM10UMZnCl1", "PM100UMZnCl2", "PM100UMZnCl1", "PM1MMZnCl2", "PM1MMZnCl1", "PMCadAGFP", "PMCadAGFP1", "PMCadAGFP01UMZnCl2", "PMCadAGFP01UMZnCl1", "PMCadAGFP1UMZnCl2", "PMCadAGFP1UMZnCl1", "PMCadAGFP10UMZnCl2", "PMCadAGFP10UMZnCl1", "PMCadAGFP100UMZnCl2", "PMCadAGFP100UMZnCl1", "PMCadAGFP1MMZnCl2", "PMCadAGFP1MMZnCl1", "PMCadAGFP2", "PMCadAGFP3", "PMCadAGFP01UMZnCl3", "PMCadAGFP01UMZnCl4", "PMCadAGFP1UMZnCl3", "PMCadAGFP1UMZnCl4", "PMCadAGFP10UMZnCl3", "PMCadAGFP10UMZnCl4", "PMCadAGFP100UMZnCl3", "PMCadAGFP100UMZnCl4", "PMCadAGFP1MMZnCl3", "PMCadAGFP1MMZnCl4", "PM4", "PM5", "PM01UMHgCl2", "PM01UMHgCl1", "PM1UMHgCl2", "PM1UMHgCl1", "PM10UMHgCl2", "PM10UMHgCl1", "PM100UMHgCl2", "PM100UMHgCl1", "PMCadAGFP4", "PMCadAGFP5", "PMCadAGFP01UMHgCl2", "PMCadAGFP01UMHgCl1", "PMCadAGFP1UMHgCl2", "PMCadAGFP1UMHgCl1", "PMCadAGFP10UMHgCl2", "PMCadAGFP10UMHgCl1", "PMCadAGFP100UMHgCl2", "PMCadAGFP100UMHgCl1"];
for i = 2:length(all_vars)
    temp_var = lower(all_vars{i});
    if contains(temp_var, 'lb')
        tempstrain = 'lb';
    elseif contains(temp_var, 'cop')
        tempstrain = 'copAGFP';
    elseif contains(temp_var, 'cad')
        tempstrain = 'cadAGFP';
    else
        tempstrain = 'PM WT';
    end
    all_strain_info_order{1, i} = tempstrain;
    if contains(temp_var, 'cu')
        temp_metal = 'cu';
    elseif contains(temp_var, 'hg')
        temp_metal = 'hg';
    elseif contains(temp_var, 'zn')
        temp_metal = 'zn';
    else
        temp_metal = 'none';
    end
    all_strain_info_order{2, i} = temp_metal;
    
    % Get the concentration with pattern matching
    conc_regexp = '\d+[um]+' ; %{1, 2}';
    conc_str = regexpi(temp_var, conc_regexp, "match");
    
    if ~isempty(conc_str)
        conc_str = conc_str{1};
        if conc_str(1) == "0"
%             disp(conc_str);
            concnum = 0.1;
        else
            concnum_regexp = '\d+';
            concnum_str = regexpi(conc_str, concnum_regexp, "match");
            concnum = str2double(concnum_str);
        end
        if contains(lower(conc_str), 'u')
            concnum = concnum/1000; % let's do everything in units of millimolar
        end
        all_strain_info_order{3, i} = concnum;
    else
        all_strain_info_order{3, i} = 0;
    end
    all_strain_info_order{4, i} = all_vars{i};
end
clear tempstrain temp_metal temp_var;
all_strain_info_order = all_strain_info_order(:, 2:end);





%% Get the fold changes plotted for each strain
strains = unique(all_strain_info_order(1, :));
metals = unique(all_strain_info_order(2, :));
concs =  unique([all_strain_info_order{3, :}]);
cop_strains = {'PM WT', 'copAGFP'};
cad_strains = {'PM WT', 'cadAGFP'};

strain_metal_combos = {'copAGFP', 'cu'; 'cadAGFP', 'zn'; 'cadAGFP', 'hg'};


% Use increasing color for more copper
colorvals = [1 0.8 0.6 0.4 0.2 0];
% Define values for patches
error_alpha = .4;
error_width_factor = .003;
figure(1); clf('reset');
figure(2); clf('reset');

% plot each combo
for combonum = 1:length(strain_metal_combos)
    temp_strain = strain_metal_combos{combonum, 1};
    temp_metal = strain_metal_combos{combonum, 2};

    figure(1);
    clf('reset');
%     subplot(1, 3, combonum);
    hold on;
    line_styles = {':', '-'};
    
    times = copAcadA_GFPtable.Timeh1(2:end);
    
    all_vecs = {};
    
    for conc_num = 1:length(concs)
        temp_conc = concs(conc_num);
        % get the inds of the columns in the table that match this
        if conc_num==1
            temp_inds = strcmpi(all_strain_info_order(1,:), temp_strain) &...
                [all_strain_info_order{3, :}]==temp_conc & ...
                strcmpi(all_strain_info_order(2, :), 'none');
        else
            temp_inds = strcmpi(all_strain_info_order(1,:), temp_strain) &...
                [all_strain_info_order{3, :}]==temp_conc & ...
                strcmpi(all_strain_info_order(2, :), temp_metal);
        end
        temp_inds = find(temp_inds);
    
        fprintf(strcat('Strain: ', temp_strain, ', metal: ', temp_metal,...
            ', conc: %g\n'), temp_conc);
        disp(length(temp_inds));
    
        temp_vecs = [];
        for ind_num = 1:length(temp_inds)
            % get the variable name
            temp_var = all_strain_info_order{4, temp_inds(ind_num)};
            temp_vecs(ind_num, :) = copAcadA_GFPtable.(temp_var)(2:end);
        end
        all_vecs{end+1} = temp_vecs;
    end
    
    % Now that we have all the vecs in order of conc, let's get the mean of
    % each and divide it by the mean of the 0 one
    fold_vecs = all_vecs;
    temp_mean_vec = mean(all_vecs{1});
    temp_mean_vec(temp_mean_vec==0) = NaN;
    for conc_num = 1:length(all_vecs)
        % get the individual vectors
        temp_vecs = all_vecs{conc_num};
        for vec_num = 1:size(temp_vecs, 1)
            temp_vecs(vec_num, :) = temp_vecs(vec_num, :)./temp_mean_vec;
        end
        fold_vecs{conc_num} = temp_vecs;
    end
    

    
    pos_ind = 5; % find(temp_mean_vec>=0, 1);
    % great, let's plot each over time
    
    for conc_num=1:length(fold_vecs)
        
        temp_conc = concs(conc_num);
        tempcolor = [colorvals(conc_num)/4, colorvals(conc_num), colorvals(conc_num)/3, 0.8];
        lname = strcat(temp_strain, " ", num2str(temp_conc), " mM ", temp_metal);

    end
    


    %%%%%%%% Let's plot the MAX fold change
    figure(2); 
    clf('reset'); 
%     subplot(1, 3, combonum);
    hold on;
    line_styles = {':', '-'};

    % If we have individual vectors divided by the 0 vec, then taking max
    % of each, then averaging
    max_fold_vals = cell(size(fold_vecs)); % these will be the max for each vector
    mean_max_fold_vals = zeros(1, length(concs));
    for conc_num = 1:length(concs)
        temp_vecs = fold_vecs{conc_num};
        for rownum = 1:size(temp_vecs, 1)
            max_fold_vals{conc_num}(rownum) = max(temp_vecs(rownum, :));
            disp(max(temp_vecs(rownum, :)));
        end
        mean_max_fold_vals(conc_num) = mean(max_fold_vals{conc_num});
    end
    % we want the first bar to be 1 since it is the one that everything
    % else is in terms of fold change of
    mean_max_fold_vals(1) = 1;

%     % Or if we have mean vector that we are dividing by mean 0 and then
%     % taking max
%     mean_max_fold_vals = zeros(1, length(concs));
%     for conc_num = 1:length(concs)
%         mean_max_fold_vals(conc_num) = max(fold_vecs{conc_num});
%     end

    % now plot each bar one at a time
    width_vals = 5*[0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1];
    b = {};
    tempconcs = concs;
    tempconcs(1) = 0.00001; % otherwise it won't plot on log scale
    for bar_num = 1:length(concs)
        b{bar_num} = bar(tempconcs(bar_num), mean_max_fold_vals(bar_num), ...
            width_vals(bar_num), ...
            'FaceColor', 'none', ...
            'EdgeColor', [0.3 0.7 0.6], ...
            'LineWidth', 2);

        temp_fold_vals = max_fold_vals{bar_num};
        % if we have individual data points
        if bar_num ~= 1
            s(i) = scatter(repelem(tempconcs(bar_num), length(temp_fold_vals)), ...
            temp_fold_vals, 150, 'filled', 'MarkerFaceColor', [0.3 0.5 0.8], ...
            'MarkerFaceAlpha', 0.7, 'XJitter', 'rand', 'XJitterWidth', 0.1);
        end
    end

    % now scatter the individual ones
%     b = bar(concs, max_fold_vals, 1.0);
%     plot(concs, max_fold_vals, 'LineWidth', 3); %, 'o', 'MarkerFaceColor', [0.2 0.8 0.4], ...
% %         'MarkerEdgeColor', 'none', 'MarkerSize', 15);
% 

    ax = gca;
    ax.XScale = 'log';
%     ax.XLim = [0.00006 1.6];
    ax.XLim = [0.000006 1.6];
%     if combonum ==3
%         ax.XLim = [0.0006 1.6];
%     end
    set(ax,'FontName','Myriad','FontSize',12,'TickDir','out','TickLength',...
    [0.02 0.02]);
    temptitle = strcat('Max Fold Change of ', temp_strain, " with ", temp_metal);
    title(temptitle,'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    xlabel(strcat(temp_metal, " Concentration (mM)"), 'FontSize', 15);
    ylabel(strcat("Max Fold Change in Fluorescence from 0 ", temp_metal),...
        'FontSize', 15);    
    % legend(plotlines, 'Interpreter', 'none', 'Location', 'best'); 
    
%     fig = gcf;
%     savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
%     % Things to do for saving
%     set(fig, 'WindowStyle', 'normal'); % in case it was docked
%     set(fig, 'Color', 'none');
%     set(gcf, 'Position', [680   385   560   420]);
%     % set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
%     drawnow;
%     pause(3);
%     savename = strcat(temp_strain, '_', temp_metal, '_max_folds_bar.svg');
%     export_fig(fig, fullfile(savefigdir, savename), '-nocrop', '-Painters');
%     
%     fig.Color = [1 1 1];
%     set(fig, 'WindowStyle', 'docked');



end

%% 01/23 Plot raw gfp values instead of max fold change

strains = unique(all_strain_info_order(1, :));
metals = unique(all_strain_info_order(2, :));
concs =  unique([all_strain_info_order{3, :}]);
cop_strains = {'PM WT', 'copAGFP'};
cad_strains = {'PM WT', 'cadAGFP'};

strain_metal_combos = {'copAGFP', 'cu'; 'cadAGFP', 'zn'; 'cadAGFP', 'hg'};


% Use increasing color for more copper
colorvals = [1 0.8 0.6 0.4 0.2 0];
% Define values for patches
error_alpha = .4;
error_width_factor = .003;
figure(1); clf('reset');
figure(2); clf('reset');

% plot each combo
for combonum = 1:length(strain_metal_combos)
    temp_strain = strain_metal_combos{combonum, 1};
    temp_metal = strain_metal_combos{combonum, 2};

    line_styles = {':', '-'};
    
    times = copAcadA_GFPtable.Timeh1(2:end);
    
    all_vecs = {};
    
    for conc_num = 1:length(concs)
        if combonum == 3 & conc_num ==6
            % no 1 mM Hg so skip
            continue;
        end
        temp_conc = concs(conc_num);
        % get the inds of the columns in the table that match this
        if conc_num==1
            temp_inds = strcmpi(all_strain_info_order(1,:), temp_strain) &...
                [all_strain_info_order{3, :}]==temp_conc & ...
                strcmpi(all_strain_info_order(2, :), 'none');
        else
            temp_inds = strcmpi(all_strain_info_order(1,:), temp_strain) &...
                [all_strain_info_order{3, :}]==temp_conc & ...
                strcmpi(all_strain_info_order(2, :), temp_metal);
        end
        temp_inds = find(temp_inds);
    
        fprintf(strcat('Strain: ', temp_strain, ', metal: ', temp_metal,...
            ', conc: %g\n'), temp_conc);
        disp(length(temp_inds));
    
        temp_vecs = [];
        for ind_num = 1:length(temp_inds)
            % get the variable name
            temp_var = all_strain_info_order{4, temp_inds(ind_num)};
            temp_vecs(ind_num, :) = copAcadA_GFPtable.(temp_var)(2:end);
        end
        if size(temp_vecs, 1) ~= 0
            all_vecs{end+1} = temp_vecs;
        end
    end
    

    % 01-23-23: Now we have all vecs in order of conc
    % Let's first try plotting all the last values?

    figure(combonum); 
    clf('reset'); 
%     subplot(1, 3, combonum);
    hold on;
    line_styles = {':', '-'};

    % If we have individual vectors divided by the 0 vec, then taking last
    % of each, then averaging
    last_vals = cell(size(all_vecs)); % these will be the max for each vector
    mean_last_vals = zeros(1, length(concs));
    for conc_num = 1:length(concs)
        if combonum == 3 & conc_num == 6
            continue;
        end
        temp_vecs = all_vecs{conc_num};
        for rownum = 1:size(temp_vecs, 1)
            last_vals{conc_num}(rownum) = temp_vecs(rownum, end);
            
        end
        mean_last_vals(conc_num) = mean(last_vals{conc_num});
    end
    

    % now plot each bar one at a time
    width_vals = 5*[0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1];
    b = {};
    tempconcs = concs;
    tempconcs(1) = 0.00001; % otherwise it won't plot on log scale
    for bar_num = 1:length(concs)
        % TODO FIX THIS PART SO IT PLOTS RIGHT 
        if combonum == 3 & bar_num == 6
            continue;
        end
        b{bar_num} = bar(tempconcs(bar_num), mean_last_vals(bar_num), ...
            width_vals(bar_num), ...
            'FaceColor', 'none', ...
            'EdgeColor', [0.3 0.7 0.6], ...
            'LineWidth', 2);

        temp_fold_vals = last_vals{bar_num};
        % if we have individual data points
        
            s(i) = scatter(repelem(tempconcs(bar_num), length(temp_fold_vals)), ...
            temp_fold_vals, 150, 'filled', 'MarkerFaceColor', [0.3 0.5 0.8], ...
            'MarkerFaceAlpha', 0.7, 'XJitter', 'rand', 'XJitterWidth', 0.1);
        
    end



    ax = gca;
    ax.XScale = 'log';
%     ax.XLim = [0.00006 1.6];
    ax.XLim = [0.000006 1.6];
    if combonum ==3
        ax.XLim = [0.000006 0.16];
    end
    set(ax,'FontName','Myriad','FontSize',12,'TickDir','out','TickLength',...
    [0.02 0.02]);
    temptitle = strcat('Fluorescence Value at end of Expt of ', temp_strain, " with ", temp_metal);
    title(temptitle,'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    xlabel(strcat(temp_metal, " Concentration (mM)"), 'FontSize', 15);
    ylabel(strcat("Experiment End Fluorescence Value ", temp_metal),...
        'FontSize', 15);    
    % legend(plotlines, 'Interpreter', 'none', 'Location', 'best'); 
    
%     fig = gcf;
%     savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
%     % Things to do for saving
%     set(fig, 'WindowStyle', 'normal'); % in case it was docked
%     set(fig, 'Color', 'none');
% %     set(gcf, 'Position', [680   385   560   420]);
%     set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
%     drawnow;
%     pause(3);
%     savename = strcat(temp_strain, '_', temp_metal, '_end_vals_bar.svg');
%     export_fig(fig, fullfile(savefigdir, savename), '-nocrop', '-Painters');
%     
%     fig.Color = [1 1 1];
%     set(fig, 'WindowStyle', 'docked');
% 
%     pause(1);

end







