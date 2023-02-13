%% Analysis and Plotting of Spot Image Analysis

%% First, get the latest analysis

spotim_dir = '/Users/anjali/Dropbox/_Shared_Patterning/Patterning_Expts/Exploratory_Expts/Fall2022';
oldFolder = cd(spotim_dir);

file_data = dir(fullfile('**', 'all_metal_expt_data_*.mat'));

if isempty(file_data)
    disp('No metal experiment data found.');
elseif length(file_data)==1
    load(fullfile(file_data(1).folder, file_data(1).name), 'allExptsTable');
elseif length(file_data)>1
    tempnames = {file_data.name};
    tempdates = NaT(1, length(tempnames));
    for i = 1:length(tempnames)
        tempdat = erase(tempnames{i}, 'all_metal_expt_data_');
        tempdat = erase(tempdat, '.mat');
        if ~isempty(tempdat)
            tempdates(i) = datetime(tempdat);
        end
    end
    % get the latest date
    [~, ind] = max(tempdates);
    % load that one
    load(fullfile(file_data(ind).folder, file_data(ind).name), 'allExptsTable');
end
clear tempdat tempnames tempdates;

%% Compile Analysis of Spot Ims & check if any unanalyzed

% Note: All this is in the allExptsTable that should be loaded

% file_data = dir(fullfile('**', '*radial_analysis.mat'));
% file_names = {file_data.name}; %cell array of file names
% fprintf('Found %d metal spot expt files \n', length(file_names));
% 
% % Now let's see if any are not yet in the table
% for i = 1:length(file_names)
%     disp(strcat("Checking for data from: ", file_names{i}));
%     mainFolder = cd(file_data(i).folder);
%     
%     % Load the major analysis
%     current = file_names{i};
%     load(current, 'polar_analysis');
%     test_name = polar_analysis(i).filename;
%     test_name = erase(test_name, 'Polar_Half_Images/');
%     % check if this image is in the allExptsTable
%     if ~any(strcmpi(test_name, allExptsTable.filename))
%         % this analysis isn't in the table yet
%         disp(strcat("Adding data from: ", file_names{i}));
%         % Load the major analysis
%         current = file_names{i};
%         load(current, 'polar_analysis');
%         % convert to table
%         t = struct2table(polar_analysis);
%         % remove extra vars
%         t2 = removevars(t, ["dist_vec", "colrad", "petrimask", "colcenter"]);
%         curr_vars = t2.Properties.VariableNames;
%         curr_height = height(allExptsTable);
%         % Add each row one by one to the table 
%         for j = 1:length(curr_vars)
%             % append all the values from the new table
%             allExptsTable.(curr_vars{j})(curr_height+1:curr_height+height(t2))...
%                 = t2.(curr_vars{j});
%         end
% % 
% 
%     end
%     % if so, clear the polar_analysis and continue
%     cd(mainFolder);
% end
% clear test_name curr_vars curr_height polar_analysis current t t2;
% cd(spotim_dir);

%% If no analysis has been done, make the table from the ground up
% We shouldn't need to do this since we loaded latest the table at the top

% % note, before doing that, we would need to turn the dist vecs to be 1000x1
% % otherwise it'll mess up the table
% for i = 1:length(file_names)
%    disp(strcat("Adding data from: ", file_names{i}));
%    mainFolder = cd(file_data(i).folder);
%    
%     % Load the major analysis
%     current = file_names{i};
%     load(current, 'polar_analysis');
%     % convert to table
%     t = struct2table(polar_analysis);
%     % remove extra vars
%     t2 = removevars(t, ["dist_vec", "colrad", "petrimask", "colcenter"]);
% 
%     % concatenate to previous table
%     if i == 1
%         allExptsTable = t2;
%     else
%         allExptsTable = [allExptsTable; t2];
%     end
% end

%% Add Promoter, Date, Metal, Conc to table
% All this is already in the table loaded at the top

%aspects to add to it: the date, the promoter, the metal, the concentration
% let's iterate over the table and add it all
% 
% for rownum = 1:height(allExptsTable)
% %     if isempty(allExptsTable.promoter{rownum})
%         tempfile = allExptsTable.filename{rownum};
%         % get the filename
%         tempfile = erase(tempfile, 'Polar_Half_Images/');
%         % get the date
%         date_pat = "202210"+ digitsPattern(2);
%         tempdate = extract(tempfile, date_pat);
%         if isempty(allExptsTable.date{rownum})
%             allExptsTable.date{rownum} = tempdate{1};
%         end
%         % get the metal
%         
%         if contains(lower(tempfile), 'cu')
%             allExptsTable.metal{rownum} = 'Cu';
%         elseif contains(lower(tempfile), 'zn')
%             allExptsTable.metal{rownum} = 'Zn';
%         else
%             allExptsTable.metal{rownum} = 'None';
%         end
%         
%     
%         % get the concentration
%     %     conc_pat = digitsPattern(2) + 'mM';
%     %     conc_pat2 = digitsPattern(2) + '_mM';
%     %     tempconc = extract(tempfile, conc_pat);
%     %     if isempty(tempconc)
%     %         tempconc = extract(tempfile, conc_pat2);
%     %     end
%     %     tempconc = erase(tempconc, 'mM');
%     %     tempconc = erase(tempconc, '_');
%     %     tempconc = str2double(tempconc);
%         
%     
%         
%         test_regexp = '\d+[_]?[m]+' ; %{1, 2}';
%         conc_str = regexpi(tempfile, test_regexp, "match");
% 
%         if isempty(conc_str)
%             allExptsTable.conc(rownum) = 0;
%             allExptsTable.metal{rownum} = 'None';
%         else
%             conc_str = erase(conc_str, 'mM');
%             conc_str = erase(conc_str, 'mm');
%             conc_str = split(conc_str, '_');
%             tempconc = str2double(conc_str(1));
%             allExptsTable.conc(rownum) = tempconc;
%         end
%         if allExptsTable.conc(rownum)==0
%             allExptsTable.metal{rownum} = 'None';
%         end
%         
%     
%         % get the gene and promoter
%     
%         if contains(lower(tempfile), 'cad')
%             allExptsTable.promoter{rownum} = 'cada';
%         elseif contains(lower(tempfile), 'cop')
%             allExptsTable.promoter{rownum} = 'copa';
%         else
%             allExptsTable.promoter{rownum} = 'None';
%         end
%     
%         if contains(lower(tempfile), 'chew')
%             allExptsTable.gene{rownum} = 'chew';
%         elseif contains(lower(tempfile), 'lrp')
%             allExptsTable.gene{rownum} = 'lrp';
%         elseif contains(lower(tempfile), 'flgm')
%             allExptsTable.gene{rownum} = 'flgm';
%         elseif contains(lower(tempfile), 'gfp')
%             allExptsTable.gene{rownum} = 'gfp';
%         elseif contains(lower(tempfile), 'wt')
%             allExptsTable.gene{rownum} = 'wt';
%         end
% 
% %     end
%     % if the row is already filled in, skip it
% end


%% Set up parameters for plotting

genes = {'gfp', 'flgm'};
concs = [0, 10, 25, 50];
temp_colors = linspecer(4);

%%
% figure(1);
% clf('reset');
% % t = tiledlayout('flow');
% for gene_num = 1:length(genes)
%     tempgene = genes{gene_num};
%     subplot(2, 1, gene_num);
%     cla('reset');
%     title(tempgene + ", sort of colony radius right v left");
%     hold on;
% %     figure(gene_num);
% %     clf('reset');
% %     t = tiledlayout(2, 2);
%     
%    
%     for conc_num = 1:length(concs)
% 
%         tempconc = concs(conc_num);
% %         nexttile;
% %         title(strcat(tempgene, num2str(tempconc)));
% %         hold on;
% 
%         tempinds = strcmpi(allExptsTable.promoter,'copa') & ...
%             allExptsTable.conc == tempconc &...
%             strcmpi(allExptsTable.gene, tempgene) &...
%             strcmpi(allExptsTable.half, 'left');
%         % get the other halves
%         tempinds_right = zeros(length(tempinds), 1);
%         for i = 1:(length(tempinds_right)-1)
%             if tempinds(i+1) == 1
%                 tempinds_right(i)=1;
%             else
%                 tempinds_right(i)=0;
%             end
%         end
%         tempinds_right = logical(tempinds_right);
% 
%         temptrajs = allExptsTable.avgd(tempinds);
%         temptrajs_right = allExptsTable.avgd(tempinds_right);
%         % plotting
%         
%         % instantiate places to store the values
%         temp_vals = zeros(length(temptrajs), 1);
% 
%         for traj_num = 1:length(temptrajs)
%             temptraj_left = temptrajs{traj_num};
%             temptraj_left = abs(gradient(movmean(temptraj_left, 10)));
%             % let's get the last place with a high gradient 
%             temptraj_ind_left = find(temptraj_left(1:500)>0.0002, 1, 'last');
% 
%             % repeat for right side
%             temptraj_right = temptrajs_right{traj_num};
%             temptraj_right = abs(gradient(movmean(temptraj_right, 10)));
%             % let's get the last place with a high gradient 
%             temptraj_ind_right = find(temptraj_right(1:500)>0.0002, 1, 'last');
% 
% %             plot(tempconc, temptraj_ind_left, 'o', "Color", temp_colors(conc_num));
% %             plot(tempconc, temptraj_ind_right, '*', "Color", temp_colors(conc_num));
% 
%             plot(tempconc, temptraj_ind_right/temptraj_ind_left, '.', "Color", temp_colors(conc_num, :), ...
%                 'MarkerSize', 12);
% %             plot(temptraj_left, "Color", [0 0.4470 0.7410]);
% %             plot(temptraj_right, '--', "Color", [0 0.2 0.3]);
% %             waitforbuttonpress;  
%             temp_vals(traj_num) = temptraj_ind_right/temptraj_ind_left;
%         end
%         % plot the mean and std
%         cv_vals = std(temp_vals)./mean(temp_vals);
%         errorbar(tempconc, mean(temp_vals), cv_vals, "-s", "MarkerSize",10,...
%             "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90]);
%         ylim([0.9, 2]);
% %         waitforbuttonpress;
%         
% 
%     end
% %     xlim([0, 600]);
% end

%% get rid of 'polar-half_images' from directory names

% for i = 1:height(allExptsTable)
%     if contains(allExptsTable.filename{i}, 'Polar_Half_Images/')
%         allExptsTable.filename{i} = erase(allExptsTable.filename{i}, ...
%             'Polar_Half_Images/');
%     end
% end

%% Write down if the image has covered plate or not

% cd(spotim_dir);
% disp('Type 1 if covering plate: ');
% for tempfilenum = 117:2:height(allExptsTable)
%     tempfile = allExptsTable.filename{tempfilenum};
%     tempfile = erase(tempfile, '_lefthalf_polarim');
%     tempfile = erase(tempfile, '_righthalf_polarim');
%     imloc = dir(fullfile('**', tempfile));
%     tempim = imread(fullfile(imloc.folder, tempfile));
%     figure(2);
%     imshow(tempim);
%     title(tempfile);
%     drawnow;
%     resp = input('Covered plate?: ');
%     if resp == 1
%         allExptsTable.covered_plate(tempfilenum) = resp;
%     else
%         allExptsTable.covered_plate(tempfilenum) = 0;
%     end
% 
% end
% 
% % fill in the left halves
% for tempfilenum = 2:2:height(allExptsTable)
%     allExptsTable.covered_plate(tempfilenum) = ...
%         allExptsTable.covered_plate(tempfilenum-1);
% end

%% Load the Mask Images from U-net and Add Analysis to the Table

% Note: This analysis should already be in table after this was run, don't
% redo when plotting

% % Iterate over each plate
% % For each, try a horizontal rolling average to fill in the gaps, then get
% % the average ring width of each ring from each column, then average over
% % the columns, to get something like [1 1.4 2] eg for 3 rings with the last
% % one being the widest
% % for the ones where it didn't cover the plate fully, get the lefthand col
% % radius versus the right one (average location of last ring edge)
% 
% cd(spotim_dir);
% for tempfilenum = 1:height(allExptsTable)
%     if isempty(allExptsTable.rng_edge_locs{tempfilenum})
%         tempfile = allExptsTable.filename{tempfilenum};
%         tempfile = strrep(tempfile, '.tif', '_skel.tif');
%         imloc = dir(fullfile('**', tempfile));
%         tempim = imread(fullfile(imloc.folder, tempfile));
%         rng_edge_vals = getRingEdges(tempim);
%         allExptsTable.rng_edge_locs{tempfilenum} = rng_edge_vals;
%     end
% end


%% get the ring widths

% for tempimnum = 1:height(allExptsTable)
%     allExptsTable.rng_widths{tempimnum} =...
%         diff(allExptsTable.rng_edge_locs{tempimnum});
% end

%% Add dist_vecs to table
% add the dist_vec to each row in the table

% cd(spotim_dir);
% file_data = dir(fullfile('**', '*radial_analysis.mat'));
% file_names = {file_data.name}; %cell array of file names
% for i = 1:length(file_names)
% 
%    disp(strcat("Loading data from: ", file_names{i}));
%    mainFolder = cd(file_data(i).folder);
%    
%     % Load the major analysis
%     current = file_names{i};
%     load(current, 'polar_analysis');
%     for j = 1:length(polar_analysis)
%         % get the image name
%         tempfile = polar_analysis(j).filename;
%         % remove polar half etc
%         tempfile = erase(tempfile, ...
%             'Polar_Half_Images/');
%         % find the index of it in allexptstable
%         table_ind = find(contains([allExptsTable.filename], tempfile));
%         % add the dist_vec to all expts table
%         allExptsTable.dist_vecs(table_ind) = {polar_analysis(j).dist_vec};
% 
%     end 
% 
% end

% clear polar_analysis

%% Now get ring edge locs and ring widths in cm
% 
% for table_ind = 1:height(allExptsTable)
%     temp_dist_vec = allExptsTable.dist_vecs{table_ind};
%     temp_rng_edge_locs = allExptsTable.rng_edge_locs{table_ind};
%     temp_rng_locs_cm = temp_rng_edge_locs;
%     for i = 1:length(temp_rng_locs_cm)
%         temp_rng_locs_cm(i) = temp_dist_vec(round(temp_rng_edge_locs(i)));
%     end
%     allExptsTable.rng_locs_cm{table_ind} = temp_rng_locs_cm;
%     allExptsTable.rng_widths_cm{table_ind} = diff(temp_rng_locs_cm);
% end
% 
% clear temp_dist_vec temp_rng_edge_locs

%% Get Colony rads in cm


% for table_ind = 1:height(allExptsTable)
%     
%     if allExptsTable.covered_plate(table_ind) ~= 1
%         temprad = max(allExptsTable.rng_locs_cm{table_ind});
%     else
%         % get the location where the plate edge is
%         temp_dist_vec = allExptsTable.dist_vecs{table_ind};
%         tempavgd = allExptsTable.avgd{table_ind};
%         temprad_ind = find(tempavgd==1, 1);
%         temprad = temp_dist_vec(temprad_ind);
%     end
%     allExptsTable.rads_cm(table_ind) = temprad;
% end
% 
% clear temp_dist_vec temp_rng_edge_locs tempavgd temprad_ind temprad






%% Get ring width values organized by day

% First, let's make a new cell with the vals for each day in the same row
% So it should go: D1, D2, D3 as rows, and across would be the 0 vals, 10mM
% vals, 25 mM vals, 50 mM vals
% we'll make one each for GFP & for flgM
genes = {'gfp', 'flgm'};
concs = [0, 10, 25, 50];

all_dates = unique(allExptsTable.date);
all_gfp_left_widths = {};
all_gfp_right_widths = {};
all_flgm_left_widths = {};
all_flgm_right_widths = {};
halves = {'left', 'right'};
for i = 1:length(all_dates)
    tempdate = all_dates{i};
    date_inds = strcmpi(allExptsTable.date,tempdate);
    % if so: we'll fill in a row for the gfp cell and the flgm cell
    curr_row = size(all_gfp_left_widths, 1)+1;
    % go through GFP & flgm
    for j = 1:length(genes)
        tempgene = genes{j};
        for k = 1:length(concs)
            tempconc = concs(k);
            for m = 1:2
                temp_half = halves{m};
                % get the inds; note if we want 0 conc, the metal will be none
                if tempconc == 0
                    tempinds = date_inds & strcmpi(allExptsTable.gene, tempgene)...
                        & strcmpi(allExptsTable.metal, 'none')...
                        & strcmpi(allExptsTable.promoter, 'copa') ...
                        & strcmpi(allExptsTable.half, temp_half)...
                        & allExptsTable.conc == tempconc;
                else
                    tempinds = date_inds & strcmpi(allExptsTable.gene, tempgene)...
                        & strcmpi(allExptsTable.metal, 'cu')...
                        & strcmpi(allExptsTable.half, temp_half)...
                        & allExptsTable.conc == tempconc;
                end
                fprintf(strcat('Date: ', tempdate, ", conc:  %d", ", gene: ", tempgene, "\n"), tempconc);
                % add the ring widths to the appropriate cell (don't overwrite previous
                temp_rng_widths = allExptsTable.rng_widths_cm(tempinds);
                if strcmpi(tempgene, 'gfp')
                    if strcmpi(temp_half, 'left')
                        all_gfp_left_widths{curr_row, k} = temp_rng_widths;
                    else
                        all_gfp_right_widths{curr_row, k} = temp_rng_widths;
                    end
                else
                    if strcmpi(temp_half, 'left')
                        all_flgm_left_widths{curr_row, k} = temp_rng_widths;
                    else
                        all_flgm_right_widths{curr_row, k} = temp_rng_widths;
                    end
                end
            end
            % continue
        end
    end
end
clear tempgene tempinds tempconc tempdate date_inds curr_row temp_rng_widths;
all_rng_width_vals = {all_gfp_left_widths, all_gfp_right_widths; all_flgm_left_widths, all_flgm_right_widths};



%% Make new matrices with the 2nd-3rd ring widths with the values normalized to day 0

% First try: Getting the 2nd-3rd ring widths, averaging each, then dividing
% each by the mean of the 0 vals for that day, then plotting

all_mean_vals = cell(2, 2);
for half_val = 1:2
    for gene_val = 1:2
        temp_all_vals = all_rng_width_vals{gene_val, half_val};
        temp_mean_vals = cell(size(temp_all_vals));
        for j = 1:size(temp_all_vals, 1)
            for k = 1:size(temp_all_vals, 2)
                temp_vals = temp_all_vals{j, k};
                for l = 1:length(temp_vals)
                    % now we are iterating over the plates
                    temp_rng_widths = temp_vals{l}(2:3);
                    temp_mean_vals{j, k}(l) = mean(temp_rng_widths);
                end
            end
        end
        all_mean_vals{gene_val, half_val} = temp_mean_vals;
    end
end

%% Plot 4l: Normalize plate by plate
% Get each plate's 2-3 ring average divided by that day's 0 plate 2-3 ring
% average
% Then scatter by concentration, color by GFP vs flgM
% can error bar at each concentration as well (the mean of all the plates)

temp_colors = {[0,144,64]/255, [224,160,0]/255};
temp_normed_vals = cell(2, 4);
for gene_num = 1:2
    % Let's try first getting only the left halves
    temp_all_vals = all_mean_vals{gene_num, 1};
    tempgene = genes{gene_num};

    
    for date_val = 1:3
        % get the 0 val for that date
        temp_0_mean_val = mean(temp_all_vals{date_val, 1});
        for conc_val = 1:4
            % Get all the plates at that conc on that day
            plate_vals = temp_all_vals{date_val, conc_val};
            % Divide each by the 0 val
            plate_vals = plate_vals/temp_0_mean_val;
            % add them to the bigger cell
            normed_vals_at_conc = temp_normed_vals{gene_num, conc_val};
            for j = 1:length(plate_vals)
                normed_vals_at_conc(end+1) = plate_vals(j);
            end
            temp_normed_vals{gene_num, conc_val} = normed_vals_at_conc;
        end
    end
    
end

% Now, we have the values of each plate's 2-3 ring widths, normalized to
% the mean of that measurement for the 0 mM CU plates of the same day
% let's scatter it and see what we get
figure(2);
clf('reset');
hold on;
for gene_num = 1:2
    disp(genes{gene_num});
    gene_mean_vals = zeros(1, 4);
    gene_cv_vals = gene_mean_vals;
    gene_std_vals = gene_mean_vals;
    tempcolor = temp_colors{gene_num};
    rng_vals_for_anova = [];
    rng_grps_for_anova = [];
    for conc_val = 1:4
        tempy = temp_normed_vals{gene_num, conc_val};
        tempx = repelem(concs(conc_val), length(tempy));
        scatter(tempx, tempy, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', tempcolor);
        gene_mean_vals(conc_val) = mean(tempy);
        gene_std_vals(conc_val) = std(tempy);
        gene_cv_vals(conc_val) = gene_std_vals(conc_val)/gene_mean_vals(conc_val);
        rng_vals_for_anova(end+1:end+length(tempy)) = tempy;
        rng_grps_for_anova(end+1:end+length(tempy)) = tempx;
        disp(length(tempx));
    end

    % Anova
    [p,~,stats] = anova1(rng_vals_for_anova, rng_grps_for_anova, 'off');
    [c,~,~,gnames] = multcompare(stats, 'Display', 'off');
    tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A") = gnames(tbl.("Group A"));
    tbl.("Group B") = gnames(tbl.("Group B"));
    disp(genes{gene_num});
    disp(tbl);
    errorbar(concs, gene_mean_vals, gene_std_vals, 'Color', tempcolor, ...
        'Marker', 'square', 'MarkerFaceColor', tempcolor, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 12);
    

end
title('Middle Ring Widths Normalized to Same-Day 0 Copper Control');

%% Format and save the plot
ax = gca;
% formatLinePlot(ax);
ax.XLim = [-5 60];
ax.XTick = [0 10 25 50];
ax.XTickLabel = {'0', '10', '25', '50'};
ax.YLim = [0.4 1.2];
ax.XLabel.String = 'Copper in Spots (mM)';
ax.YLabel.String = 'Normalized Ring Width';
% ax.YTick = [1 2 3 4];
ticklengthval = [0 0];
tickdirval = 'out';
keeplinewidths=true;
axfontsize = 24; smallticks = false; 
titlesize = 12;
plotgrid = false;
legend('', '', '', '', 'GFP', '', '', '', '', 'flgM');
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

% 
% fig = gcf;
% savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
% % Things to do for saving
% set(fig, 'WindowStyle', 'normal'); % in case it was docked
% set(fig, 'Color', 'none');
% % set(gcf, 'Position', [680   385   560   420]);
% set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
% drawnow;
% export_fig(fig, fullfile(savefigdir, 'Fig4_m.svg'), '-nocrop', '-Painters');
% 
% fig.Color = [1 1 1];
% set(fig, 'WindowStyle', 'docked');

%% Plot 16c right panel: Get the ring width values normalized to GFP instead of 0

temp_normed_vals = cell(2, 4);
% Let's try first getting only the left halves
temp_flgm_vals = all_mean_vals{2, 1};
temp_gfp_vals = all_mean_vals{1, 1};
tempgene = 'flgM';

for gene_num = 1:2
    temp_all_vals = all_mean_vals{gene_num, 1};
    tempgene = genes{gene_num};
    for date_val = 1:3
        
        for conc_val = 1:4
            % get the gfp val for that date & conc
            temp_gfp_mean_val = mean(temp_gfp_vals{date_val, conc_val});
            % Get all the plates at that conc on that day
            plate_vals = temp_all_vals{date_val, conc_val};
            % Divide each by the 0 val
            plate_vals = plate_vals/temp_gfp_mean_val;
            % add them to the bigger cell
            normed_vals_at_conc = temp_normed_vals{gene_num, conc_val};
            for j = 1:length(plate_vals)
                normed_vals_at_conc(end+1) = plate_vals(j);
            end
            temp_normed_vals{gene_num, conc_val} = normed_vals_at_conc;
        end
    end
end

figure(3);
clf('reset');
hold on;
for gene_num = 1:2
    gene_mean_vals = zeros(1, 4);
    gene_cv_vals = gene_mean_vals;
    gene_std_vals = gene_mean_vals;
    tempcolor = temp_colors{gene_num};
    rng_vals_for_anova = [];
    rng_grps_for_anova = [];
    for conc_val = 1:4
        tempy = temp_normed_vals{gene_num, conc_val};
        tempx = repelem(concs(conc_val), length(tempy));
        scatter(tempx, tempy, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', tempcolor);
        gene_mean_vals(conc_val) = mean(tempy);
        gene_std_vals(conc_val) = std(tempy);
        gene_cv_vals(conc_val) = gene_std_vals(conc_val)/gene_mean_vals(conc_val);
        rng_vals_for_anova(end+1:end+length(tempy)) = tempy;
        rng_grps_for_anova(end+1:end+length(tempy)) = tempx;
    
    end
    
    % Anova
    [p,~,stats] = anova1(rng_vals_for_anova, rng_grps_for_anova, 'off');
    [c,~,~,gnames] = multcompare(stats, 'Display', 'off');
    tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A") = gnames(tbl.("Group A"));
    tbl.("Group B") = gnames(tbl.("Group B"));
    disp(genes{gene_num});
    disp(tbl);
    errorbar(concs, gene_mean_vals, gene_std_vals, 'Color', tempcolor, ...
        'Marker', 'square', 'MarkerFaceColor', tempcolor, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 12);

end

title('Middle Ring Widths Normalized to Same-Day GFP Control');

%% Now let's format & save it
ax = gca;
% formatLinePlot(ax);
ax.XLim = [-5 60];
ax.XTick = [0 10 25 50];
ax.XTickLabel = {'0', '10', '25', '50'};
ax.YLim = [0.3 1.1];
ax.XLabel.String = 'Copper in Spots (mM)';
ax.YLabel.String = 'Normalized Ring Width';
% ax.YTick = [1 2 3 4];
ticklengthval = [0 0];
tickdirval = 'out';
keeplinewidths=true;
axfontsize = 24; smallticks = false; 
titlesize = 12;
plotgrid = false;
legend('', '', '', '', 'GFP', '', '', '', '', 'flgM', ...
    'location', 'best');
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);


% fig = gcf;
% savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
% % Things to do for saving
% set(fig, 'WindowStyle', 'normal'); % in case it was docked
% set(fig, 'Color', 'none');
% % set(gcf, 'Position', [680   385   560   420]);
% set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
% drawnow;
% export_fig(fig, fullfile(savefigdir, 'FigS16_c_coppernormgfp_plot.svg'), '-nocrop', '-Painters');
% 
% fig.Color = [1 1 1];
% set(fig, 'WindowStyle', 'docked');

%% Plot 16c left panel: Plot the values NOT normalized to any controls

temp_normed_vals = cell(2, 4);
% Let's try first getting only the left halves
temp_flgm_vals = all_mean_vals{2, 1};
temp_gfp_vals = all_mean_vals{1, 1};
tempgene = 'flgM';

for gene_num = 1:2
    temp_all_vals = all_mean_vals{gene_num, 1};
    tempgene = genes{gene_num};
    for date_val = 1:3
        
        for conc_val = 1:4
            
            % Get all the plates at that conc on that day
            plate_vals = temp_all_vals{date_val, conc_val};
            
            % add them to the bigger cell
            normed_vals_at_conc = temp_normed_vals{gene_num, conc_val};
            for j = 1:length(plate_vals)
                normed_vals_at_conc(end+1) = plate_vals(j);
            end
            temp_normed_vals{gene_num, conc_val} = normed_vals_at_conc;
        end
    end
end

figure(4);
clf('reset');
hold on;
for gene_num = 1:2
    gene_mean_vals = zeros(1, 4);
    gene_cv_vals = gene_mean_vals;
    gene_std_vals = gene_mean_vals;
    tempcolor = temp_colors{gene_num};
    rng_vals_for_anova = [];
    rng_grps_for_anova = [];
    for conc_val = 1:4
        tempy = temp_normed_vals{gene_num, conc_val};
        tempx = repelem(concs(conc_val), length(tempy));
        scatter(tempx, tempy, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', tempcolor);
        gene_mean_vals(conc_val) = mean(tempy);
        gene_std_vals(conc_val) = std(tempy);
        gene_cv_vals(conc_val) = gene_std_vals(conc_val)/gene_mean_vals(conc_val);
        rng_vals_for_anova(end+1:end+length(tempy)) = tempy;
        rng_grps_for_anova(end+1:end+length(tempy)) = tempx;
    
    end
    
    % Anova
    [p,~,stats] = anova1(rng_vals_for_anova, rng_grps_for_anova, 'off');
    [c,~,~,gnames] = multcompare(stats, 'Display', 'off');
    tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A") = gnames(tbl.("Group A"));
    tbl.("Group B") = gnames(tbl.("Group B"));
    disp(genes{gene_num});
    disp(tbl);
    errorbar(concs, gene_mean_vals, gene_std_vals, 'Color', tempcolor, ...
        'Marker', 'square', 'MarkerFaceColor', tempcolor, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 12);

end

title('Middle Ring Widths Not Normalized to Controls');

%% Now let's format & save it
ax = gca;
% formatLinePlot(ax);
ax.XLim = [-5 60];
ax.XTick = [0 10 25 50];
ax.XTickLabel = {'0', '10', '25', '50'};
% ax.YLim = [0.3 1.1];
ax.XLabel.String = 'Copper in Spots (mM)';
ax.YLabel.String = 'Mean Middle Ring Width (pixels)';
% ax.YTick = [1 2 3 4];
ticklengthval = [0 0];
tickdirval = 'out';
keeplinewidths=true;
axfontsize = 24; smallticks = false; 
titlesize = 12;
plotgrid = false;
legend('', '', '', '', 'GFP', '', '', '', '', 'flgM', ...
    'location', 'best');
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);


fig = gcf;
savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
% set(gcf, 'Position', [680   385   560   420]);
set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
drawnow;
export_fig(fig, fullfile(savefigdir, 'FigS16_c_notnorm_plot.svg'), '-nocrop', '-Painters');

fig.Color = [1 1 1];
set(fig, 'WindowStyle', 'docked');



%% Get radius values organized by day

% First, let's make a new cell with the vals for each day in the same row
% So it should go: D1, D2, D3 as rows, and across would be the 0 vals, 10mM
% vals, 25 mM vals, 50 mM vals
% we'll make one each for GFP & for flgM
genes = {'gfp', 'flgm'};
concs = [0, 10, 25, 50];

all_dates = unique(allExptsTable.date);
all_gfp_left_rads = {};
all_gfp_right_rads = {};
all_flgm_left_rads = {};
all_flgm_right_rads = {};
halves = {'left', 'right'};
for i = 1:length(all_dates)
    tempdate = all_dates{i};
    date_inds = strcmpi(allExptsTable.date,tempdate);
    % if so: we'll fill in a row for the gfp cell and the flgm cell
    curr_row = size(all_gfp_left_rads, 1)+1;
    % go through GFP & flgm
    for j = 1:length(genes)
        tempgene = genes{j};
        for k = 1:length(concs)
            tempconc = concs(k);
            for m = 1:2
                temp_half = halves{m};
                % get the inds; note if we want 0 conc, the metal will be none
                if tempconc == 0
                    tempinds = date_inds & strcmpi(allExptsTable.gene, tempgene)...
                        & strcmpi(allExptsTable.metal, 'none')...
                        & strcmpi(allExptsTable.promoter, 'copa') ...
                        & strcmpi(allExptsTable.half, temp_half)...
                        & allExptsTable.conc == tempconc;
                else
                    tempinds = date_inds & strcmpi(allExptsTable.gene, tempgene)...
                        & strcmpi(allExptsTable.metal, 'cu')...
                        & strcmpi(allExptsTable.half, temp_half)...
                        & allExptsTable.conc == tempconc;
                end
                fprintf(strcat('Date: ', tempdate, ", conc:  %d", ", gene: ", tempgene, "\n"), tempconc);
                disp(length(find(tempinds)));
                % add the ring widths to the appropriate cell (don't overwrite previous
                temp_rng_rads = allExptsTable.rads_cm(tempinds);


                % Store the radii

                if strcmpi(tempgene, 'gfp')
                    if strcmpi(temp_half, 'left')
                        all_gfp_left_rads{curr_row, k} = temp_rng_rads;
                    else
                        all_gfp_right_rads{curr_row, k} = temp_rng_rads;
                    end
                else
                    if strcmpi(temp_half, 'left')
                        all_flgm_left_rads{curr_row, k} = temp_rng_rads;
                    else
                        all_flgm_right_rads{curr_row, k} = temp_rng_rads;
                    end
                end
            end
            % continue
        end
    end
end
% clear tempgene tempinds tempconc tempdate date_inds curr_row temp_rng_locs temp_rng_rads;
all_rad_vals = {all_gfp_left_rads, all_gfp_right_rads; all_flgm_left_rads, all_flgm_right_rads};

%% Make new matrices with the colony 'radius' with the values normalized to day 0

% First try: Getting the rads, dividing
% each by the mean of the 0 vals for that day, then plotting

all_mean_vals = cell(2, 2);
for half_val = 1:2
    % left vs right half
    for gene_val = 1:2
        % now iterating over the two genes
        temp_all_vals = all_rad_vals{gene_val, half_val};
        temp_mean_vals = cell(size(temp_all_vals));
        for j = 1:size(temp_all_vals, 1)
            % now iterating over the expts (rows)
            for k = 1:size(temp_all_vals, 2)
                % now iterating over the columns dates)
                temp_vals = temp_all_vals{j, k};
                for l = 1:length(temp_vals)
                    % now we are iterating over the and getthing the radius
                    % of each
                    temp_mean_vals{j, k}(l) = temp_vals(l);
                end
            end
        end
        all_mean_vals{gene_val, half_val} = temp_mean_vals;
    end
end





%% Plot 16b - Plotting the rads plate by plate NOT normalized
% Get each plate's rad divided by that day's 0 copper rad avg
% Then scatter by concentration, color by GFP vs flgM
% can error bar at each concentration as well (the mean of all the plates)

temp_normed_vals = cell(2, 4);
for gene_num = 1:2
    % Let's try first getting only the left halves
    temp_all_vals = all_mean_vals{gene_num, 1};
    tempgene = genes{gene_num};

    
    for date_val = 1:3
        % get the 0 val for that date
        temp_0_mean_val = mean(temp_all_vals{date_val, 1});
        for conc_val = 1:4
            % Get all the plates at that conc on that day
            plate_vals = temp_all_vals{date_val, conc_val};
            % DON't normalize the plate vals
            % add them to the bigger cell
            normed_vals_at_conc = temp_normed_vals{gene_num, conc_val};
            for j = 1:length(plate_vals)
                normed_vals_at_conc(end+1) = plate_vals(j);
            end
            temp_normed_vals{gene_num, conc_val} = normed_vals_at_conc;
        end
    end
    
end

% Now, we have the values of each plate's 2-3 ring widths, normalized to
% the mean of that measurement for the 0 mM CU plates of the same day
% let's scatter it and see what we get
figure(3);
clf('reset');
hold on;
for gene_num = 1:2
    gene_mean_vals = zeros(1, 4);
    gene_cv_vals = gene_mean_vals;
    gene_std_vals = gene_mean_vals;
    tempcolor = temp_colors{gene_num};
    rng_vals_for_anova = [];
    rng_grps_for_anova = [];
    for conc_val = 1:4
        tempy = temp_normed_vals{gene_num, conc_val};
        tempx = repelem(concs(conc_val), length(tempy));
        scatter(tempx, tempy, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', tempcolor);
        gene_mean_vals(conc_val) = mean(tempy);
        gene_std_vals(conc_val) = std(tempy);
        gene_cv_vals(conc_val) = gene_std_vals(conc_val)/gene_mean_vals(conc_val);
        rng_vals_for_anova(end+1:end+length(tempy)) = tempy;
        rng_grps_for_anova(end+1:end+length(tempy)) = tempx;

    end

    % Anova
    disp('Anova for Normed Rad Vals');
    [p,~,stats] = anova1(rng_vals_for_anova, rng_grps_for_anova, 'off');
    [c,~,~,gnames] = multcompare(stats, 'Display', 'off');
    tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A") = gnames(tbl.("Group A"));
    tbl.("Group B") = gnames(tbl.("Group B"));
    disp(genes{gene_num});
    disp(tbl);
    errorbar(concs, gene_mean_vals, gene_std_vals, 'Color', tempcolor, ...
        'Marker', 'square', 'MarkerFaceColor', tempcolor, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 12);
    

end
title('Colony Radii');

%% Now let's format & save it
ax = gca;
% formatLinePlot(ax);
ax.XLim = [-5 60];
ax.XTick = [0 10 25 50];
ax.XTickLabel = {'0', '10', '25', '50'};
ax.YLim = [0 5];
ax.XLabel.String = 'Copper in Spots (mM)';
ax.YLabel.String = 'Colony Radius (cm)';
% ax.YTick = [1 2 3 4];
ticklengthval = [0 0];
tickdirval = 'out';
keeplinewidths=true;
axfontsize = 24; smallticks = false; 
titlesize = 12;
plotgrid = false;
legend('', '', '', '', 'GFP', '', '', '', '', 'flgM');
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);


fig = gcf;
savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
% set(gcf, 'Position', [680   385   560   420]);
set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
drawnow;
export_fig(fig, fullfile(savefigdir, 'Fig_S15_copRads_normedto0.svg'), '-nocrop', '-Painters');

fig.Color = [1 1 1];
set(fig, 'WindowStyle', 'docked');


%% Save updated table
cd(spotim_dir);
newfilename = strcat('all_metal_expt_data_', date, '.mat');
save(newfilename, 'allExptsTable', '-v7.3');
cd(oldFolder);
%% OG save

% disp('Completed compiling timelapse data. Saving to current folder...');
% % SAVE THE STRUCT TO A FILE
% cd(spotim_dir);
% % Save the measurements in a .mat file in corresponding folder, use a
% % version compatible with large variables. Note: future updates should
% % be saved with the date in the name.
% filename = 'all_metal_expt_data.mat';
% save(filename, 'allExptsTable', '-v7.3');
% cd(oldFolder);
% disp('Done.');


%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
function rng_edge_vals = getRingEdges(tempim)
    % in this function, from a ring mask output by unet, get approximate
    % ring edge locations

    % get the average number of ring edges
    mean_ring_edge_num = round(mean(sum(tempim, 1)));
    edge_vals = zeros(round(mean_ring_edge_num)+5, 1000);
    %iterate over the columns
    for i = 1:size(tempim, 2)
        tempcol = tempim(:, i);
        tempinds = find(tempcol);
        edge_vals(1:length(tempinds), i) = tempinds;
    end
    edge_vals(edge_vals==0) = NaN;
    plot(nanmean(edge_vals, 2));
    rng_edge_vals = nanmean(edge_vals, 2);
    rng_edge_vals = rng_edge_vals(~isnan(rng_edge_vals));
    % remove ones that are too close
    while any(diff(rng_edge_vals)<15)
        if any(diff(rng_edge_vals)<15)
            % likely not to be true rings at 15 pixel distance
            repeat_ring_ind = find(diff(rng_edge_vals)<15);
    
            mean_edge_loc = mean([rng_edge_vals(repeat_ring_ind), ...
                rng_edge_vals(repeat_ring_ind+1)]);
            rng_edge_vals(repeat_ring_ind:repeat_ring_ind+1) = mean_edge_loc;
            rng_edge_vals = unique(rng_edge_vals);
        end
    end
    % if the last one is at a lower ind than the first one, remove 
    end_ind = length(rng_edge_vals);
    if rng_edge_vals(end_ind)<rng_edge_vals(end_ind-1)
        rng_edge_vals = rng_edge_vals(end_ind-1);
    end
end

function formatLinePlot(ax_all, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
%     ax.FontName = 'Myriad';

    if ~exist('plotgrid', 'var')
        plotgrid = false;
    end
    
    
    for axind = 1:length(ax_all)
        ax = ax_all(axind);
        if ~strcmpi(get(ax_all, 'Type'), 'ConfusionMatrixChart') &...
                ~strcmpi(get(ax_all, 'Type'), 'heatmap')

            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.XColor = [0.2 0.2 0.2]; 
            ax.YColor = [0.2 0.2 0.2];
            ax.TickDir = tickdirval;    
            ax.TickLength = ticklengthval;

            if smallticks
                ax.FontSize = 14;
                if ~isempty(ax.XLabel)
                    ax.XLabel.FontSize = 24;
                end
                if ~isempty(ax.YLabel)
                    ax.YLabel.FontSize = 24;
                end
            else
                % Change both x & y labels, and tick mark fonts
                ax.FontSize = axfontsize;
            end

            % get lines
            if ~keeplinewidths
                if length(ax.Children) == 1
                    linetype = ax.Children.Type;
                    if strcmpi(linetype, 'errorbar')
                        errorlines = ax.Children;
                        errorlines.LineWidth = 3;
                        errorlines.CapSize = 1;
                    elseif strcmpi(linetype, 'line')
                        lines = ax.Children;
                        lines.LineWidth = 3;
                    end
                elseif length(ax.Children)> 1
                    for j = 1:length(ax.Children)
                        temp_part = ax.Children(j);
                        temptype = ax.Children(j).Type;
                        if strcmpi(temptype, 'errorbar')
                            temp_part.LineWidth = 3;
                            temp_part.CapSize = 1;
                        elseif strcmpi(temptype, 'line')
                            temp_part.LineWidth = 3;
                        elseif strcmpi(temptype, 'scatter')
                            % for now don't do anything

                        end
                    end
                end
            end


            % Check if legend, if so, format
            if ~isempty(ax.Legend)
                ax.Legend.Color = [1 1 1];
                ax.Legend.EdgeColor = [0.9 0.9 0.9];
                ax.Legend.FontSize = 12;
                ax.Legend.LineWidth = 2;
                ax.Legend.TextColor = ax.XColor;
            end

            if ~isempty(ax.Title)
                ax.Title.Color = ax.XColor;
                ax.Title.FontSize = titlesize;
                ax.Title.FontWeight = 'bold';
            end

            % These options can be changed if we want a background color/grid
            if plotgrid
                ax.Color = [0.9098    0.9098    0.9098];
                ax.XGrid = 'on'; 
                ax.YGrid = 'on';
                ax.GridColor = [0.95 0.95 0.95];
                ax.GridAlpha = 1;
                ax.MinorGridColor = [0.95 0.95 0.95];
                ax.MinorGridLineStyle = '-';
                ax.MinorGridAlpha = 1;
            end
        end
    end
end