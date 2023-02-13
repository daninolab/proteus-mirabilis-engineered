%% Combo Sensor AUC 05-24-21

% Make table from combo sensor data, fit model to classify into the 9
% classes, & then calculate multi-class AUC on it
%% Load latest Data
cd(comboExptDir);
allExptData = loadLatestComboSensorData();

%% Making AUC table--Make the Full Table

fnames = fieldnames(allExptData);

% % Disp fnames with inds if desired
% for i = 1:length(fnames)
%     disp(strcat(num2str(i), " ", fnames{i})); 
% end

% Decide fields to keep
fnames = fnames([3, 4, 13, 14, 23, 24, 15, 17:18, 22, 26:27, 29:33]);
for i = 1:length(allExptData)
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

%% Setting up class bins thing--adding a column to tdata
iptg_bounds = [0, 1, 5];
ara_bounds = [0, 0.1, 0.2];
class_bins = [1, 2, 3; 4, 5, 6; 7, 8, 9];
for i = 1:height(tdata)
    tempiptg = tdata.iptgs(i);
    tempara = tdata.aras(i);
    iptg_ind = find(tempiptg>=iptg_bounds, 1, 'last');
    ara_ind = find(tempara>=ara_bounds, 1, 'last');
    tdata.classbins(i) = class_bins(iptg_ind, ara_ind);
end
    
%% Set Up Params to use & classes
%%%% Loop over param combos & store AUC

% set up the param combos to use
all_param_combs = {};
last_var_val = length(tdata.Properties.VariableNames)-1;
for i = 1:12
    temp_combs = nchoosek([5:last_var_val], i);
    for j = 1:size(temp_combs, 1)
        all_param_combs{end+1} = temp_combs(j, :);
    end
end


classbins = [1, 2, 3, 4, 5, 6, 7, 8, 9];
classLabels = categorical({'0-1 mM IPTG, 0-0.09% ara',...
    '1-4.9 mM IPTG, 0-0.09% ara', '5-10 mM IPTG, 0-0.09% ara', ...
    '0-1 mM IPTG, 0.1-0.19% ara', '1-4.9 mM IPTG, 0.1-0.19% ara',...
    '5-10 mM IPTG, 0.1-0.19% ara', '0-1 mM IPTG, 0.2% ara', ...
    '1-4.9 mM IPTG, 0.2% ara', '5-10 mM IPTG, 0.2% ara'});


% set up places to store aucs
all_aucs = zeros(1, length(all_param_combs)); 

% store warnings
all_warns = all_aucs;
y_ord = ordinal(tdata.classbins);

%% Calculate AUCs

% loop over the param combos
for j = 1:length(all_param_combs)
    if rem(j, 50)==0
        fprintf('Combo %g of %g\n', j, length(all_param_combs));
    end
    vars_to_use = all_param_combs{j}; %variables in tdata2
    % Get the x data for these variables
    x = tdata{:, vars_to_use};
    
    % Reset the warning so we can catch it later:
    warning('');
    
    % Fit multinomial logistic model
    b_ord = mnrfit(x, y_ord, 'model', 'ordinal');
    pihat_ord = mnrval(b_ord,x,'model','ordinal');
    
    % Calculate multiClassAUC & store
    pClass = pihat_ord;
    classLabel = double(y_ord);
    all_aucs(j) = multiClassAUC(pClass, classLabel);
    
    % Check if we had a warning
    if ~isempty(lastwarn)
        all_warns(i, j, k) = 1;
    end
end
    
%% Get Nice Variable Combo Names
varnames = tdata.Properties.VariableNames;
varcombo_names = cell(1, length(all_param_combs));
for j = 1:length(varcombo_names)
    tempstr = strjoin(varnames(all_param_combs{j}));
    tempstr = strrep(tempstr, '_', " ");
    varcombo_names{j} = tempstr;
end


    

%% Get dataset size
dataset_size = zeros(3, 3);
for i = 1:9
    dataset_size(i) = length(find(tdata.classbins==i));
end
dataset_size = dataset_size';
disp(dataset_size);
    
    

%% Plot confusion matrix & save things to a struct

comboGeneAUCData = struct('gene', {}, 'numsamples', {}, ...
    'numsamples_10mm', {}, 'b_ord', {}, 'vars_used', {}, ...
    'varnames_used', {}, 'weights', {},...
    'est_cats', {}, 'y_ord', {},  'max_auc', {});

% Calculate the best-fit model & fit
% Then use calculations to make CM

% Get the max AUC where the warning is none--get index of parameters &
% class binning
temp_aucs = all_aucs;
tempwarns = all_warns;
temp_aucs(tempwarns==1) = NaN;
[max_auc, max_combo_ind] = max(temp_aucs(:, :));


comboGeneAUCData(1).max_auc = max_auc;
comboGeneAUCData.vars_used = all_param_combs{max_combo_ind};
comboGeneAUCData.varnames_used = varnames(all_param_combs{max_combo_ind});

% Get the parameters for that to use from tdata2
vars_to_use = all_param_combs{max_combo_ind};

% Get the x data for these variables
x = tdata{:, vars_to_use};

class_bins = class_bins;
classLabels = classLabels;


% Use these class labels to sort targets
% y_ord = ordinal(tdata.classbins, {'1', '2', '3', '4', '5', '6', '7', '8', '9'}, [], class_bins);
y_ord = ordinal(tdata.classbins);
comboGeneAUCData.y_ord = y_ord;

% Fit multinomial logistic model
b_ord = mnrfit(x, y_ord, 'model', 'ordinal');
pihat_ord = mnrval(b_ord,x,'model','ordinal');
comboGeneAUCData.b_ord = b_ord;


% Calculate multiClassAUC & store
pClass = pihat_ord;
classLabel = double(y_ord);
tempauc = multiClassAUC(pClass, classLabel);


[~, est_cats] = max(pihat_ord, [], 2);
est_cats = ordinal(est_cats);
comboGeneAUCData.est_cats = est_cats;

% Get weights
temp_weights = zeros(1, length(vars_to_use));
for j = 1:length(vars_to_use)
    weight_val = b_ord(length(b_ord)-length(vars_to_use)+j);
    temp_weights(j) = weight_val;
end
comboGeneAUCData.weights = temp_weights;



new_classlabels = cell(1, length(classLabels));
for i = 1:length(new_classlabels)
    new_classlabels{i} = string(classLabels(i));
end
figure(2);
clf('reset');
% plotconfusion(y, est_cats);
% cm = confusionchart(y_ord, est_cats);
cm = confusionchart(double(y_ord), double(est_cats));
% titletext = strcat("Combo Sensor Prediction Accuracy binning into ", ...
%     strjoin(string(classLabels), ", "));
titletext = {"Combo Sensor Prediction Accuracy binning into ", ...
    new_classlabels{:}, ", "};

title({titletext{:}, sprintf("mAUC: %g", tempauc)});
ax = gca;
% ax.FontSize = 18;


%% Save AUC

disp('Saving AUCs...');
filename = strcat('combo_AUCs_', date, '.mat');
save(filename, 'comboGeneAUCData', '-v7.3');
disp('Done.');
    
%% FUNCTIONS

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
        
        if isdatetime(tempdata)
            % it's the date
            tempdate = string(tempdata);
            tempdate = strrep(tempdate, '-00', '-20');
            for k = 1:numims
                tempvals{k, j} = tempdate;
            end
        end
    end

end
