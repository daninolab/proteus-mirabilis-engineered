%% Goal: Calculate & Plot AUCs for Fitting Curves for Single Gene Data

% Load Single Gene Data
cd(singleExptDir);
allExptData = loadLatestSingleGeneData();




    
%% Making AUC table--Make the Full Table

fnames = fieldnames(allExptData);
% Decide fields to keep


% TRYING ON 9/11/21: Removing certain fields which are duplicates etc
% to do: edit this & keep the mean intens, sector cv, min intens, sector
% radius in CM only!, sector area, sector mean local CV, and make sure all
% are bin style 1 only before we make this latex table

fnames = fnames([13, 14, 4, 27:29, 31, 35, 42], :);

gene_list = {'wt', 'gfp', 'umod', 'lrp', 'chew', 'flgm', 'flia'};
num_aug = length(allExptData(1).aug_col_areas{1});
for i = 1:length(allExptData)
    if isempty(allExptData(i).scanfolder)
        continue;
    end
    tempvals = mkExptCellAug(allExptData, i, fnames, num_aug);
    % now that we got all the values we want into a cell, let's make it a table
    temptbl = cell2table(tempvals, 'VariableNames', fnames);
    if i == 1
        tdata = temptbl;
    else
        % concatenate
        tdata = [tdata; temptbl];
    end
end
% remove outlier rows
inds = tdata.outlier==0;
tdata = tdata(inds, :);
% remove outlier column
tdata = removevars(tdata, 'outlier');
%% Multinomial VUCs
%%%% Loop over genes, class organizations, param combos & store AUC

% set up the param combos to use
all_param_combs = {};
for i = 1:9
    temp_combs = nchoosek([3:7], i);
    for j = 1:size(temp_combs, 1)
        all_param_combs{end+1} = temp_combs(j, :);
    end
end

% set up the two class organizations desired
classbins = {[0, 0.9, 5, 10], [0, 0.09, 0.9, 10]};
classLabels1 = categorical({'0-0.9 mM', '1-5 mM', '5-10 mM'});
classLabels2 = categorical({'0-0.09 mM', '0.1-0.9 mM', '1-10 mM'});
classLabelsAll = {classLabels1, classLabels2};
% set up places to store aucs
all_aucs = zeros(7, length(all_param_combs), 2); 

% store warnings
all_warns = all_aucs;
% rows = genes, columns = parameter combination, z direction = which of the
% class bins were used

% loop over genes
for i = 1:length(gene_list)
    gene_curr = gene_list{i};
    % Obtain all the values for that gene
    inds = (strcmpi(gene_curr, tdata.genes));
    fprintf('Gene %g has %g samples\n', i, length(find(inds)));
    tdata2 = tdata(inds, [1, 2, 4:8]); %, 7, 8]);

    % loop over the param combos
    for j = 1:length(all_param_combs)
        if rem(j, 50)==0
            fprintf('Gene %g combo %g of %g\n', i, j, length(all_param_combs));
        end
        vars_to_use = all_param_combs{j}; %variables in tdata2
        % Get the x data for these variables
        x = tdata2{:, vars_to_use};
        
        % loop over the two class organization
        for k = 1:2
            
            % Reset the warning so we can catch it later:
            warning('');
            
            class_bins = classbins{k};
            classLabels = classLabelsAll{k};
            
            % Use these class labels to sort targets
            y_ord = ordinal(tdata2.iptgs, {'1', '2', '3'}, [], class_bins);

            % Fit multinomial logistic model
            b_ord = mnrfit(x, y_ord, 'model', 'ordinal');
            pihat_ord = mnrval(b_ord,x,'model','ordinal');

            % Calculate multiClassAUC & store
            pClass = pihat_ord;
            classLabel = double(y_ord);
            all_aucs(i, j, k) = multiClassAUC(pClass, classLabel);
            
            % Check if we had a warning
            if ~isempty(lastwarn)
                all_warns(i, j, k) = 1;
            end
            
        end
    end
end

%% Get Nice Variable Combo Names
varnames = tdata2.Properties.VariableNames;
varcombo_names = cell(1, length(all_param_combs));
for j = 1:length(varcombo_names)
    tempstr = strjoin(varnames(all_param_combs{j}));
%     tempstr = replace(tempstr, "percent_area", "% Area");
%     tempstr = replace(tempstr, "colrads_cm", "Radii");
%     tempstr = replace(tempstr, "avgd", "Mean Intens");
%     tempstr = replace(tempstr, "cvs", "Mean CV");
%     tempstr = replace(tempstr, "stdevs", "Mean Stdev");
    varcombo_names{j} = tempstr;
end



%% Get Max Aucs
% get max gene_aucs
temp_aucs = all_aucs;
temp_aucs(all_warns==1) = NaN;

for i = 1:7
    bin_ind = 1;
    
    [maxauc, maxind] = nanmax(temp_aucs(i, :, 1));
    [maxauc2, maxind2] = nanmax(temp_aucs(i, :, 2));
    
%     [maxauc, maxind] = max(all_aucs(i, :, 1));
%     [maxauc2, maxind2] = max(all_aucs(i, :, 2));
    if maxauc2>maxauc
        maxauc = maxauc2;
        maxind = maxind2;
        bin_ind = 2;
    end
    disp(gene_list(i));
    disp(varcombo_names{maxind});
    disp(maxauc);
    fprintf('Max AUC Bin Class: %g\n', bin_ind);
end

%% Print max aucs nicely
aucs_cell = cell(7, 4);
for i = 1:7
    bin_ind = 1;
    [maxauc, maxind] = max(all_aucs(i, :, 1));
    [maxauc2, maxind2] = max(all_aucs(i, :, 2));
%     if maxauc2>maxauc
%         maxauc = maxauc2;
%         maxind = maxind2;
%         bin_ind = 2;
%     end
    bin_ind = 1;
    
    aucs_cell{i, 1} = gene_list{i};
    aucs_cell{i, 2} = varcombo_names{maxind};
    aucs_cell{i, 3} = maxauc;
    aucs_cell{i, 4} = bin_ind;
end
    

%% Count # samples per gene
for i = 1:length(gene_list)
    gene_curr = gene_list{i};
%     inds = (tdata.iptgs==0 | tdata.iptgs==1 | tdata.iptgs == 0.1 | tdata.iptgs == 10) & ...
%         (strcmpi(gene_curr, tdata.genes));
    inds = (strcmpi(gene_curr, tdata.genes));
    fprintf(strcat("Gene ", gene_curr, " has %g samples \n"),length(find(inds)));
    inds = (tdata.iptgs==10 & ...
        (strcmpi(gene_curr, tdata.genes)));
    fprintf(strcat("Gene ", gene_curr, " has %g 10mm plates \n"),length(find(inds)));
end

%% Plot the best fit ones
i = 7;
gene_curr = gene_list{i};
% Obtain all the values for that gene
inds = (strcmpi(gene_curr, tdata.genes));
tdata2 = tdata(inds, [1, 2, 4:9]); %, 7, 8]);

% Get the max AUC where the warning is none--get index of parameters &
% class binning
temp_aucs = all_aucs(i, :, :);
tempwarns = all_warns(i, :, :);
temp_aucs(tempwarns==1) = NaN;
[max_auc, maxind] = max(temp_aucs(:, :, 1));
[max_auc2, maxind2] = max(temp_aucs(:, :, 2));
if max_auc>max_auc2
    max_combo_ind = maxind;
    max_bin_ind = 1;
else
    max_combo_ind = maxind2;
    max_bin_ind = 2;
end
% Get the parameters for that to use from tdata2
vars_to_use = all_param_combs{max_combo_ind};
% Get the x data for these variables
x = tdata2{:, vars_to_use};

class_bins = classbins{max_bin_ind};
classLabels = classLabelsAll{max_bin_ind};


% Use these class labels to sort targets
y_ord = ordinal(tdata2.iptgs, {'1', '2', '3'}, [], class_bins);

% Fit multinomial logistic model
b_ord = mnrfit(x, y_ord, 'model', 'ordinal');
pihat_ord = mnrval(b_ord,x,'model','ordinal');

% Calculate multiClassAUC & store
pClass = pihat_ord;
classLabel = double(y_ord);
tempauc = multiClassAUC(pClass, classLabel);

% Plot Confus Matreix
[~, est_cats] = max(pihat_ord, [], 2);
est_cats = ordinal(est_cats);
figure(2);
clf('reset');
% plotconfusion(y, est_cats);
cm = confusionchart(double(y_ord), double(est_cats));
titletext = strcat(gene_curr, " Prediction Accuracy binning into ", strjoin(string(classLabels), ", "));

title({titletext, sprintf("mAUC: %g", tempauc), gene_curr});
ax = gca;
ax.FontSize = 18;
% ax.Title.FontSize = 14;

%% Get the weights of the best fit ones

allvarnames = fnames(5:end);
all_weights = cell(8, length(allvarnames)+1);
for i = 1:length(allvarnames)
    all_weights{1, i+1} = allvarnames{i};
end
all_p = all_weights;

for i = 1:7
    gene_curr = gene_list{i};
    all_weights{i+1, 1} = gene_curr;
    all_p{i+1, 1} = gene_curr;
    % Obtain all the values for that gene
    inds = (strcmpi(gene_curr, tdata.genes));
    tdata2 = tdata(inds, [1, 2, 4:8]);

    % Get the max AUC where the warning is none--get index of parameters &
    % class binning
    temp_aucs = all_aucs(i, :, :);
    tempwarns = all_warns(i, :, :);
    temp_aucs(tempwarns==1) = NaN;
    [max_auc, maxind] = max(temp_aucs(:, :, 1));
    

    
    % For the sake of figure, using bin style 1 for everything
    max_combo_ind = maxind;
    max_bin_ind = 1;
    
    % Get the parameters for that to use from tdata2
    vars_to_use = all_param_combs{max_combo_ind};
    var_names_used = varnames(all_param_combs{max_combo_ind});
    
    % Get the x dat & class bins a for these variables
    x = tdata2{:, vars_to_use};
    class_bins = classbins{max_bin_ind};
    classLabels = classLabelsAll{max_bin_ind};

    % Use these class labels to sort targets
    y_ord = ordinal(tdata2.iptgs, {'1', '2', '3'}, [], class_bins);

    % Fit multinomial logistic model
    [b_ord, dev, stats] = mnrfit(x, y_ord, 'model', 'ordinal');
    
    % Use b_ord to get weights & store
    for j = 1:length(var_names_used)
        weight_val = b_ord(length(b_ord)-length(var_names_used)+j);
        var_ind = find(strcmpi(var_names_used{j}, allvarnames));
        all_weights{i+1, var_ind+1} = weight_val;
        all_p{i+1, var_ind+1} = stats.p(length(b_ord)-length(var_names_used)+j);
    end
    
    disp(length(b_ord));
    fprintf("P value is %g\n", [stats.p]);
    
end

%% %% Get the weights of the best fit ones AND SAVE IN A STRUCT

singleGeneAUCData = struct('gene', {}, 'numsamples', {}, ...
    'numsamples_10mm', {}, 'b_ord', {}, 'vars_used', {}, 'weights', {},...
    'est_cats', {}, 'y_ord', {},  'max_auc', {});
for i = 1:length(gene_list)
    singleGeneAUCData(i).gene = gene_list{i};
end

% Get num Samples per gene
for i = 1:length(gene_list)
    gene_curr = gene_list{i};
%     inds = (tdata.iptgs==0 | tdata.iptgs==1 | tdata.iptgs == 0.1 | tdata.iptgs == 10) & ...
%         (strcmpi(gene_curr, tdata.genes));
    inds = (strcmpi(gene_curr, tdata.genes));
    num_samples = length(find(inds));
    singleGeneAUCData(i).numsamples = num_samples;
    inds = (tdata.iptgs==10 & ...
        (strcmpi(gene_curr, tdata.genes)));
    num_samples = length(find(inds));
    singleGeneAUCData(i).numsamples_10mm = num_samples;
end


allvarnames = fnames(5:end);
all_weights = cell(8, length(allvarnames)+1);
for i = 1:length(allvarnames)
    all_weights{1, i+1} = allvarnames{i};
end

% Iterate over the genes; find and store values in struct as we do so
for i = 1:7
    % skip umod for now
    if i == 3
        continue;
    end
    
    gene_curr = gene_list{i};
    all_weights{i+1, 1} = gene_curr;
    % Obtain all the values for that gene
    inds = (strcmpi(gene_curr, tdata.genes));
    tdata2 = tdata(inds, [1, 2, 4:9]);

    % Get the max AUC where the warning is none--get index of parameters &
    % class binning
    temp_aucs = all_aucs(i, :, :);
    tempwarns = all_warns(i, :, :);
    temp_aucs(tempwarns==1) = NaN;
    [max_auc, maxind] = max(temp_aucs(:, :, 1));
    % STORE IN THE STRUCT
    singleGeneAUCData(i).max_auc = max_auc;
    
    % For the sake of figure, using bin style 1 for everything
    max_combo_ind = maxind;
    max_bin_ind = 1;
    
    % Get the parameters for that to use from tdata2
    vars_to_use = all_param_combs{max_combo_ind};
    var_names_used = varnames(all_param_combs{max_combo_ind});
    singleGeneAUCData(i).vars_used = var_names_used;
    
    % Get the x dat & class bins a for these variables
    x = tdata2{:, vars_to_use};
    class_bins = classbins{max_bin_ind};
    classLabels = classLabelsAll{max_bin_ind};

    % Use these class labels to sort targets
    y_ord = ordinal(tdata2.iptgs, {'1', '2', '3'}, [], class_bins);
    singleGeneAUCData(i).y_ord = y_ord;
    
    % Fit multinomial logistic model
    b_ord = mnrfit(x, y_ord, 'model', 'ordinal');
    singleGeneAUCData(i).b_ord = b_ord;
    
    % Use b_ord to get weights & store
    temp_weights = zeros(1, length(var_names_used));
    for j = 1:length(var_names_used)
        weight_val = b_ord(length(b_ord)-length(var_names_used)+j);
        var_ind = find(strcmpi(var_names_used{j}, allvarnames));
        all_weights{i+1, var_ind+1} = weight_val;
        temp_weights(j) = weight_val;
    end
    singleGeneAUCData(i).weights = temp_weights;
    
    % Get est cats
    pihat_ord = mnrval(b_ord,x,'model','ordinal');

    % Calculate multiClassAUC & store
    pClass = pihat_ord;
    classLabel = double(y_ord);
    tempauc = multiClassAUC(pClass, classLabel);

    % Get the estimated/predicted categories
    [~, est_cats] = max(pihat_ord, [], 2);
%     est_cats = ordinal(est_cats);
    
    singleGeneAUCData(i).est_cats = est_cats;
    
    disp(length(b_ord));
    
end


%% convert confusion matrix to percentages
class_vals = cm.NormalizedValues;
% convert to percent of each row
total_ims = sum(sum(class_vals));
total_per_class = sum(class_vals, 2);
class_percents = class_vals/total_ims;
% 
% for i = 1:length(total_per_class)
%     class_percents(i, :) = class_percents(i, :)/total_per_class(i);
% end

% convert to integers
class_perc_ints = round(100*class_percents);
% class_perc_ints(class_perc_ints==0) = 1;
figure(3);
clf('reset');
classLabels = categorical({'0-0.9 mM', '1-5 mM', '5-10 mM'});
cm2 = confusionchart(int8(class_perc_ints), classLabels);
titletext = strcat(gene_curr, " Prediction Accuracy binning into ", strjoin(string(classLabels), ", "));
title(titletext);

%% umoD AUCs specifically--unaugmented--using inocfeatures
inocFtPath = '/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/Single_Gene_Analyses/inoc_features_structs/umoD_inoc_features_06-03-21.mat';
load(inocFtPath, 'inocFeatures');
% Make the main table
fnames = fieldnames(allExptData);
% % Decide fields to keep
% To only use inoc msmts:
fnames = fnames([13, 14, 4, 27, 41]);
% To use aug msmts from earlier:
% fnames = fnames([13, 14, 4, 27:35, 42], :);

gene_list = {'wt', 'gfp', 'umod', 'lrp', 'chew', 'flgm', 'flia'};
for i = 1:length(allExptData)
    if isempty(allExptData(i).scanfolder)
        continue;
    end
    tempvals = mkExptCell(allExptData, i, fnames);
    % now that we got all the values we want into a cell, let's make it a table
    temptbl = cell2table(tempvals, 'VariableNames', fnames);
    if i == 1
        tdata = temptbl;
    else
        % concatenate
        tdata = [tdata; temptbl];
    end
end
% remove outlier rows
inds = tdata.outlier==0;
tdata = tdata(inds, :);
% remove outlier column
tdata = removevars(tdata, 'outlier');

% Keep only the umod rows
tdata_umod_inds = (strcmpi('umod', tdata.genes));
tdata_umod = tdata(tdata_umod_inds, :); %, 7, 8]);

umod_fnames = fieldnames(inocFeatures);
disp(umod_fnames);
umod_fnames = umod_fnames([6:9, 18]);
% Add the fields from Marian struct
for i = 1:length(inocFeatures)
    disp(i);
   for j = 1:length(inocFeatures(i).ImagesList)
       % Load the image name
       imname = inocFeatures(i).ImagesList{j};
       imname = erase(imname, '_polarim');
       if contains(imname, '.jpg')
           imname = erase(imname, '.tif');
       end
       
       % If the imname contains the date, get rid of it
       if strcmpi(inocFeatures(i).ExptFolder(1:3), imname(1:3))
           if ~contains(imname, '10-23-20')
               imname = imname(10:end); 
           else
               imname = strrep(imname, '10-23-20_', '10_23_20_');
           end
       end
       
       
       % Find the row in tdata
%        imind = find(contains(tdata_umod.imfiles, imname));
       imind = find(strcmpi(imname, tdata_umod.imfiles));
       
       % Append each of the following fields: max intensity, min intensity,
       % mean intensity, inoc area, inoc minor axis length
       if ~isempty(imind)
           for k = 1:length(umod_fnames)
               tdata_umod.(umod_fnames{k})(imind) = inocFeatures(i).(umod_fnames{k}){j};
           end
       
       end
    
    
   end
end

tdata_umod = tdata_umod([1:7, 9:end], :); 
disp('done setting up for fitting umod');

%%
% Fit the model
% set up the param combos to use
% Subtract 6: 'properties', row, variables, genes, iptgs, & imfiles
num_params = length(fieldnames(tdata_umod))-6; 

% nchoosek 4:8 if using inoc msmts only
% 4:12 if using all the aug msmts
all_param_combs = {};
for i = 1:num_params
    temp_combs = nchoosek([5:9], i); 
    for j = 1:size(temp_combs, 1)
        all_param_combs{end+1} = temp_combs(j, :);
    end
end

% set up the two class organizations desired
classbins = {[0, 0.9, 5, 10], [0, 0.09, 0.9, 10]};
classLabels1 = categorical({'0-0.9 mM', '1-5 mM', '5-10 mM'});
classLabels2 = categorical({'0-0.09 mM', '0.1-0.9 mM', '1-10 mM'});
classLabelsAll = {classLabels1, classLabels2};
% set up places to store aucs
all_aucs = zeros(length(all_param_combs), 2); 

% store warnings
all_warns = all_aucs;
% rows = genes, columns = parameter combination, z direction = which of the
% class bins were used

% loop over the param combos
for j = 1:length(all_param_combs)
    if rem(j, 50)==0
        fprintf('Gene %g combo %g of %g\n', i, j, length(all_param_combs));
    end
    vars_to_use = all_param_combs{j}; %variables in tdata2
    % Get the x data for these variables
    x = tdata_umod{:, vars_to_use};
    if iscell(x)
        x = cell2mat(x);
    end

    % loop over the two class organization
    for k = 1:2

        % Reset the warning so we can catch it later:
        warning('');

        class_bins = classbins{k};
        classLabels = classLabelsAll{k};

        % Use these class labels to sort targets
        y_ord = ordinal(tdata_umod.iptgs, {'1', '2', '3'}, [], class_bins);

        % Fit multinomial logistic model
        b_ord = mnrfit(x, y_ord, 'model', 'ordinal');
        pihat_ord = mnrval(b_ord,x,'model','ordinal');

        % Calculate multiClassAUC & store
        pClass = pihat_ord;
        classLabel = double(y_ord);
        all_aucs(j, k) = multiClassAUC(pClass, classLabel);
        if isempty(multiClassAUC(pClass, classLabel))
            disp('found');
%             break;
        end
        % Check if we had a warning
        if ~isempty(lastwarn)
            all_warns(j, k) = 1;
        end

    end
end

%

% Get Nice Variable Combo Names
varnames = tdata_umod.Properties.VariableNames;
varcombo_names = cell(1, length(all_param_combs));
for j = 1:length(varcombo_names)
    tempstr = strjoin(varnames(all_param_combs{j}));
    varcombo_names{j} = tempstr;
end


%% Plot umoD Conf Matrix & SAVE IN STRUCT

gene_curr = 'umoD';
% Obtain all the values for that gene

tdata2 = tdata_umod;

% Get the max AUC where the warning is none--get index of parameters &
% class binning
temp_aucs = all_aucs;
tempwarns = all_warns;
temp_aucs(tempwarns==1) = NaN;
[max_auc, maxind] = max(temp_aucs(:, 1));
%%
if max_auc==1
    temp_auc_bin1 = temp_aucs(:, 1);
    max_auc = max(temp_auc_bin1(temp_auc_bin1 ~= 1));
    maxind = find(temp_auc_bin1==max_auc, 1, 'last');
end
% [max_auc2, maxind2] = max(temp_aucs(:, 2));
% if max_auc>max_auc2
%     max_combo_ind = maxind;
%     max_bin_ind = 1;
% else
%     max_combo_ind = maxind2;
%     max_bin_ind = 2;
% end

% USING BIN STYLE 1 ONLY
max_combo_ind = maxind;
max_bin_ind = 1;

singleGeneAUCData(3).max_auc = max_auc;


% Get the parameters for that to use from tdata2
vars_to_use = all_param_combs{max_combo_ind};
var_names_used = varnames(all_param_combs{max_combo_ind});
singleGeneAUCData(3).vars_used = var_names_used;

% Get the x data for these variables
x = tdata2{:, vars_to_use};
if iscell(x)
    x = cell2mat(x);
end

class_bins = classbins{max_bin_ind};
classLabels = classLabelsAll{max_bin_ind};


% Use these class labels to sort targets
y_ord = ordinal(tdata2.iptgs, {'1', '2', '3'}, [], class_bins);
singleGeneAUCData(3).y_ord = y_ord;

% Fit multinomial logistic model
[b_ord, dev, stats] = mnrfit(x, y_ord, 'model', 'ordinal');
pihat_ord = mnrval(b_ord,x,'model','ordinal');

singleGeneAUCData(3).b_ord = b_ord;

% Calculate multiClassAUC & store
pClass = pihat_ord;
classLabel = double(y_ord);
tempauc = multiClassAUC(pClass, classLabel);

% Plot Confuse Matrix
[~, est_cats] = max(pihat_ord, [], 2);
est_cats = ordinal(est_cats);
singleGeneAUCData(3).est_cats = est_cats;

var_names_used = varnames(all_param_combs{max_combo_ind});
% Get the weights
% Use b_ord to get weights & store
temp_weights = zeros(1, length(var_names_used));
all_p = temp_weights;
for j = 1:length(var_names_used)
    weight_val = b_ord(length(b_ord)-length(var_names_used)+j);
    var_ind = find(strcmpi(var_names_used{j}, allvarnames));
    temp_weights(j) = weight_val;
    all_p(j) = stats.p(length(b_ord)-length(var_names_used)+j);
end
singleGeneAUCData(3).weights = temp_weights;


figure(2);
clf('reset');
% plotconfusion(y, est_cats);
cm = confusionchart(double(y_ord), double(est_cats));
titletext = strcat(gene_curr, " Prediction Accuracy binning into ", strjoin(string(classLabels), ", "));

title({titletext, sprintf("mAUC: %g", tempauc), gene_curr});
ax = gca;
ax.FontSize = 18;
% ax.Title.FontSize = 14;

%% Save the single gene auc

disp('Saving AUCs...');
filename = strcat('single_gene_AUCs_', date, '.mat');
save(filename, 'singleGeneAUCData', '-v7.3');
disp('Done.');












%% Functions

function tempvals = mkExptCellAug(allExptData, i, fnames, num_aug)
    numims = length(allExptData(i).imfiles);
    tempvals = cell(num_aug*numims, length(fnames));
    
    for j = 1:length(fnames)
        % Get the data for that row, fieldname
        tempdata = allExptData(i).(fnames{j});
        if j < 5 %| j > 12 
            tempdata = repelem(tempdata, num_aug);
%         elseif j > 12
%             tempdata = allExptData(i).(fnames{j});
%             tempdata = repelem(tempdata, num_aug);
        else
            % if j > 4, now we need to like, break out tempdata into individual
            tempdata = (cell2mat(tempdata))';
        end
        
        % Get the length
        if ~isstr(tempdata) & length(tempdata)>1
        % If the length is > 1, iterate over to populate the cell
            for k = 1:(num_aug*numims)
                if ~iscell(tempdata)
                    tempvals{k, j} = tempdata(k);
                else
                    tempvals{k, j} = tempdata{k};
                end

            end
        end
    end

end

% NON AUGMENTED

function tempvals = mkExptCell(allExptData, i, fnames)
    numims = length(allExptData(i).imfiles);
    tempvals = cell(numims, length(fnames));
    
    for j = 1:length(fnames)
        tempdata = allExptData(i).(fnames{j});
        % Get the length
        if ~isstr(tempdata) & length(tempdata)>1
        % If the length is > 1, iterate over to populate the cell
            for k = 1:(numims)
                if ~iscell(tempdata)
                    tempvals{k, j} = tempdata(k);
                else
                    tempvals{k, j} = tempdata{k};
                end

            end
        end
    end
end