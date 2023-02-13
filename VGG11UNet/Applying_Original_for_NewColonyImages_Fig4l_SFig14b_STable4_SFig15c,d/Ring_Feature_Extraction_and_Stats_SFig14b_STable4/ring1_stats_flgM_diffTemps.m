%% Statistical testing of flgM ring #1 measurments 
% (as extracted using postprocessed VGG-11 U-Net predicted masks)
% at 0 vs. 10 mM iptg for 37C, 36C, and 34C.
% Statistics used in paper include: 
    % mean width
    % width standard deviation 
    % p-values from 2-tailed t-tests
        % (significance of difference in mean width at 0 vs. 10 mM iptg,
        % at each temperature)
    % power estimation of above 2-tailed t-tests

%%
clear all;
clc;

% get struct files that have ring msmts 
% for each flgM img within the temperature-varying expts
temperatureFolder = '/Users/marianshaw/Dropbox/U-Net_Revision_Expts/temp_sensitvity_power_analysis';
condition_msmtsFolder = strcat(temperatureFolder,'/msmt_structs_per_condition');
struct_list = dir(fullfile(condition_msmtsFolder,'*.mat'));
struct_names = {struct_list.name}';
num_structs = length(struct_names);

% set up new struct for storing statistics of select ring #
new_struct = struct('Original_Msmts_File',[],...
                    'MeanWidths_Pix',[],...
                    'MeanWidths_CM',[],...
                    'OverallMeanWidth_CM',[],...
                    'OverallSTDWidth_CM',[],...
                    'OverallSEMWidth_CM',[]);

% need the polar analysis files (from the original plattening process)
% for converting measurements from pixel to cm
polarFolder = strcat(temperatureFolder,'/polar_analysis_files');
pol_list = dir(fullfile(polarFolder, '*.mat'));
pol_names = {pol_list.name}';
all_img_names = [];
all_dist_vecs = [];
for p = 1:length(pol_names)
    file_p = pol_names{p};
    path_p = strcat(polarFolder, '/', file_p);
    load(path_p); % saved under variable "polar_analysis"
    img_list = {polar_analysis.filename}';
    dist_vec = {polar_analysis.dist_vec}';
    all_img_names = [all_img_names; img_list];
    all_dist_vecs = [all_dist_vecs; dist_vec];
end

% set the ring # we want to look at 
% eg '2' is the 1st ring (after inoculum)
ring_idx = 2; 

% get the statistics of this ring across all images 
% at a given iptg & temperature
for n = 1:num_structs
    
    % grab 1 struct 
    this_struct = struct_names{n};
    struct_path = strcat(condition_msmtsFolder, '/', this_struct);
    load(struct_path); % this will be saved under the variable 'temperature_struct' 
    
    % figure out the temperature & iptg 
    split_name = split(this_struct,'_');
    iptg_str = split_name{2};
    iptg = str2double(erase(iptg_str,'i'));
    temp = split_name{3};
    
    % get columns of imgNames & mean ring width info 
    ring_widths = {temperature_struct.MeanDistancesBtwnBoundaries}';
    img_names = {temperature_struct.ImageName}';
    num_imgs = length(ring_widths);
    
    % in a vector, store the mean widths of the selected ring for all imgs  
    % at this iptg & temp 
    cm_vec = []; 
    pix_vex = [];
    for i = 1:num_imgs
        % determine image name & expt date 
        % in order to get the corresponding pixel-->cm info
        img_i = img_names{i};
        img_i_name = erase(img_i,'_polarim');
        p_idx = find(strcmp(all_img_names,img_i_name));
        dist_vec_p = all_dist_vecs{p_idx};
        % CONVERT TO CM: 
        % get the mean difference between rows 
        meandiff = mean(diff(dist_vec_p));
        % Multiply the width in pixels by meandiff to get approximate width in cm
        width_pix = ring_widths{i}(ring_idx);
        pix_vex = [pix_vex width_pix];
        width_cm = width_pix*meandiff;
        cm_vec = [cm_vec width_cm];
    end

    % stats - calculate for the msmts already converted to cm
    mean_width = mean(cm_vec);
    std_width = std(cm_vec);
    sem_width = std_width/sqrt(num_imgs);
    
    % save all info in new struct
    new_struct(n).Original_Msmts_File = {this_struct};
    new_struct(n).MeanWidths_Pix = {pix_vex};
    new_struct(n).MeanWidths_CM = {cm_vec};
    new_struct(n).OverallMeanWidth_CM = mean_width;
    new_struct(n).OverallSTDWidth_CM = std_width;
    new_struct(n).OverallSEMWidth_CM = sem_width;
end

% save the stats struct
stats_struct_path = strcat(temperatureFolder, '/diffTemps_flgM_ring1_stats.mat');
save(stats_struct_path,'new_struct');

%% SIGNIFICANCE 
% these p-values are indicated w/ astericks on Fig. S14b
% and exact values are specified in Fig. S14b caption

% make vectors for the measurments [cm] at all conditions
all_msmts_cm = {new_struct.MeanWidths_CM};
cm_34C_0i = cell2mat(all_msmts_cm{1});
cm_34C_10i = cell2mat(all_msmts_cm{4});
cm_36C_0i = cell2mat(all_msmts_cm{2});
cm_36C_10i = cell2mat(all_msmts_cm{5});
cm_37C_0i = cell2mat(all_msmts_cm{3});
cm_37C_10i = cell2mat(all_msmts_cm{6});

% determine variance type of 0 vs. 10 mM iptg ring msmts @ each temp
[h1,p1,ci1,stats1] = vartest2(cm_37C_10i,cm_37C_0i); % h = 1 (unequal var)
[h2,p2,ci2,stats2] = vartest2(cm_36C_10i,cm_36C_0i); % h = 1 (unequal var)
[h3,p3,ci3,stats3] = vartest2(cm_34C_10i,cm_34C_0i); % h = 0 (equal var)

% 2-tailed t-tests
[h_37,p_37,ci_37,stats_37] = ttest2(cm_37C_0i,cm_37C_10i,'Vartype','unequal'); % P = 2.4216e-06
[h_36,p_36,ci_36,stats_36] = ttest2(cm_36C_0i,cm_36C_10i,'Vartype','unequal'); % P = 0.00359
[h_34,p_34,ci_34,stats_34] = ttest2(cm_34C_0i,cm_34C_10i,'Vartype','equal'); % P = 0.001027

%% POWER ANALYSIS 
% extpowerStudent is a MATLAB fileexchange
% downloaded from:
% https://www.mathworks.com/matlabcentral/fileexchange/23184-extpowerstudent

% these power estimations are tabulated in Supplementary Table 4. 

num_tails = 2;
sig_level = 0.05;
pwr_37 = extpowerStudent(stats_37.tstat,stats_37.df,num_tails,sig_level); % pwr = 1.0000
pwr_36 = extpowerStudent(stats_36.tstat,stats_36.df,num_tails,sig_level); % pwr = 0.9982
pwr_34 = extpowerStudent(stats_34.tstat,stats_34.df,num_tails,sig_level); % pwr = 0.9917

