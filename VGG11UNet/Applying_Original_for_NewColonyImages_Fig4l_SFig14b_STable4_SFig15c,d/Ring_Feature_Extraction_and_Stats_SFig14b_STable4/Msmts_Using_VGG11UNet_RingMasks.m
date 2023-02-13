%% Msmts_Using_VGG11UNet_RingMasks
% Full script for taking ring-related measurements using 
% ring boundary masks (the already skeletonized versions) 
% predicted the VGG-11 U-Net model 

%% Some initial settings
clear all;
clc;

% Set temperature, strain, & iptg concentration 
set_temp = 34; % one of: 34, 36, 37
set_iptg = 10; % one of: 0, 10
set_strain = 'flgM'; 

% Set folders of images & masks
temperatureFolder = '/Users/marianshaw/Dropbox/U-Net_Revision_Expts/temp_sensitvity_power_analysis';
condition_msmtsFolder = strcat(temperatureFolder,'/msmt_structs_per_condition');
imgsFolder =  strcat(temperatureFolder,'/',int2str(set_temp),'C','/',int2str(set_iptg),'i','/',set_strain,'_',int2str(set_iptg),'i_',int2str(set_temp),'C');
uNetFolder = strcat(temperatureFolder,'/All_U-Net_RingBoundaries');
predMasksFolder = strcat(uNetFolder,'/predictions');
skelMasksFolder = strcat(uNetFolder,'/skel_predictions');

% List of images to analyze here
img_list = dir(fullfile(imgsFolder, '*_polarim.tif'));
img_names = {img_list.name}';
num_imgs = length(img_list);

% Set up struct for storing measurements
valueHolder = {};
[temperature_struct] = create_temperature_struct(valueHolder);

%% Make sub-folders to store:
% (1) red versions of predicted masks (not skeletonized)
red_preds_sub = strcat(imgsFolder, '/uNet_preds_red');
if ~isfolder(red_preds_sub)
    mkdir(red_preds_sub);
end
% (2) red predicted masks overlayed
red_preds_overs_sub = strcat(imgsFolder, '/uNet_preds_red_over');
if ~isfolder(red_preds_overs_sub)
    mkdir(red_preds_overs_sub);
end
% (3) post-processed masks
postproc_sub = strcat(imgsFolder, '/postproc');
if ~isfolder(postproc_sub)
    mkdir(postproc_sub);
end
% (4) red versions of post-processed masks
red_postproc_sub = strcat(imgsFolder, '/postproc_red');
if ~isfolder(red_postproc_sub)
    mkdir(red_postproc_sub);
end
% (5) red post-processed masks + yellow ignored columns overlayed 
red_postproc_over_sub = strcat(imgsFolder, '/postproc_red_over');
if ~isfolder(red_postproc_over_sub)
    mkdir(red_postproc_over_sub);
end

%% Iterate through the images:
% refine the masks, 
% and take measurments

for n = 1:num_imgs
    
    img_name = img_names{n};
    
    % fill struct with basic file info
    temperature_struct(n).ImageName = img_name; 
    temperature_struct(n).Strain = set_strain; 
    temperature_struct(n).Temperature = set_temp; 
    temperature_struct(n).IPTG = set_iptg; 
    
    % read in image & masks
    [this_img, pred_mask, skel_mask] = read_img_and_masks(imgsFolder, predMasksFolder, skelMasksFolder, img_name);
    
    % determine which columns have the num_boundaries_mode
    [num_boundary_pixels, labeled_mask, num_boundaries_mode, cols_to_consider, num_cols_consider] = cols_with_boundaries_mode(skel_mask);
    
    % get rid of cols that dont have num_boundaries_mode
    [reduced_mask, vertical_count] = reduce_mask_for_calcs(num_boundaries_mode, num_cols_consider, cols_to_consider, labeled_mask);
    
    % insert initial counts into msmt struct
    [temperature_struct] = insert_initial_counts(temperature_struct, n, num_boundary_pixels, num_boundaries_mode, vertical_count);
    
    % determine the locations of ring boundaries in each column
    [locations] = get_boundary_locations(num_boundaries_mode, num_cols_consider, cols_to_consider, reduced_mask);
    
    % remove outlier boundary pixels (i.e. in btwn true rings)
    [cleaned_mask, final_cols_to_consider, final_num_cols_consider] = remove_outlier_pixels(reduced_mask, num_boundaries_mode, locations, cols_to_consider);
    
    % determine the locations of the final ring boundaries in each column
    [final_locations] = get_boundary_locations(num_boundaries_mode, final_num_cols_consider, final_cols_to_consider, cleaned_mask);
    
    % FEATURE EXTRACTION USING THE FINAL CLEANED MASK
    % calculate distance btwn boundaries
    [overall_min_max_mean_dist, means_stds_per_dist] = distance_btwn_boundaries(num_boundaries_mode, final_num_cols_consider, final_cols_to_consider, final_locations);
    % calculate how boundary locations change
    [std_boundary_locations] = calc_boundary_slopes(final_locations, num_boundaries_mode);
    
    % fill struct w/ feature msmts
    [temperature_struct] = insert_temperature_msmts(n, temperature_struct, final_num_cols_consider, overall_min_max_mean_dist, means_stds_per_dist, std_boundary_locations);

    % create red versions of masks & overlays
    [redPred, redFin, redPredOver, redFinOver_Yel] = red_masks_and_overlays(pred_mask, cleaned_mask, this_img, final_cols_to_consider);
    
    % save masks & overlays
    red_pred_path = strcat(red_preds_sub, '/', img_name);
    imwrite(redPred, red_pred_path, 'tif');
    red_pred_over_path = strcat(red_preds_overs_sub, '/', img_name);
    imwrite(redPredOver, red_pred_over_path, 'tif');
    postproc_path = strcat(postproc_sub, '/', img_name);
    imwrite(cleaned_mask, postproc_path, 'tif');
    red_postproc_path = strcat(red_postproc_sub, '/', img_name);
    imwrite(redFin, red_postproc_path, 'tif');
    red_postproc_over_path = strcat(red_postproc_over_sub, '/', img_name);
    imwrite(redFinOver_Yel, red_postproc_over_path, 'tif');
end


% save final struct
dateFormat = 'mm-dd-yy';
timeFormat = 'HH-MM-SS';
currentDate = datestr(now,dateFormat);
currentTime = datestr(now,timeFormat);
currentDateTime = strcat(currentDate, '_', currentTime);
newStructName = strcat(set_strain,'_',int2str(set_iptg),'i_',int2str(set_temp),'C_rbStruct_',currentDateTime, '.mat');
newStructPath = strcat(condition_msmtsFolder, '/', newStructName);
save(newStructPath,'temperature_struct');
