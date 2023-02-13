% Make fig 4 overlays of pred & ground truth flgM masks

% skel_imdir = '/Users/anjali/Downloads/vgg11_unet_flgm_with_bibi_results_111021/';
% skel_imdir = '/Users/anjali/Downloads/new_unet_preds_122021/';
skel_imdir = '/Users/anjali/Dropbox/_Shared_Patterning/_final_swarm_code_NCB/VGG11UNet/Reimplementation_for_ConditionSwitchingBoundaries_Figs_4j_S11e/Fig_4j_masks/';
skel_im_data = dir(fullfile(skel_imdir, '*.tif'));
skel_imfiles = {skel_im_data.name};
polarimfiles = skel_imfiles;
for i = 1:length(polarimfiles)
    polarimfiles{i} = strrep(skel_imfiles{i}, '_skel_ep59_postproc.tif', '.tif');
end

polarimdir = '/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/Single_Gene_Expts/datasets_for_col_segmentation/flgm_subset/test_ims/testims_masked/uncropped_boundary_masks/';
% OLD
% output_imdir = '/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels';

% NEW
output_imdir = '/Users/anjali/Dropbox/Shared_Swarming_Manuscript/NatChemBio_Reviews_Round2/Updated_Figs_Captions/New_Panels/';

gtfiles = skel_imfiles;
for i = 1:length(skel_imfiles)
    gtfiles{i} = strrep(polarimfiles{i}, '.tif', '_testim_boundarymask_uncrop.tif');
end
% Iterate over the images, load the masks, load the image, make overlay

for i = 1:length(skel_imfiles)
    tempimname = polarimfiles{i};
    temp_skelim = im2double(imread(fullfile(skel_imdir, skel_imfiles{i})));
    temp_gtim = im2double(imread(fullfile(polarimdir, gtfiles{i})));

    
    tempim = im2double(imread(fullfile(polarimdir, polarimfiles{i})));
    temp_skelmask = logical(temp_skelim);
    temp_gtmask = logical(temp_gtim);
    % From the ground truth ims, remove small objects
    temp_gtmask = bwareaopen(temp_gtmask,20);
    
    % Can't see anything: dilate
    se = strel('disk', 3, 8);
    temp_skelmask = imdilate(temp_skelmask, se);
    temp_gtmask = imdilate(temp_gtmask, se);
    
    % Make overlay images
    gt_maskim_r = tempim; 
    gt_maskim_g = gt_maskim_r;
    gt_maskim_b = gt_maskim_g;
    gt_maskim_r(temp_gtmask) = 0.1;
    gt_maskim_g(temp_gtmask) = 0.5;
    gt_maskim_b(temp_gtmask) = 0.8;
    gt_maskim = cat(3, gt_maskim_r, gt_maskim_g, gt_maskim_b);
    
    pred_maskim_r = tempim; 
    pred_maskim_g = pred_maskim_r;
    pred_maskim_b = pred_maskim_g;
    pred_maskim_r(temp_skelmask) = 0.8;
    pred_maskim_g(temp_skelmask) = 0.5;
    pred_maskim_b(temp_skelmask) = 0.1;
    pred_maskim = cat(3, pred_maskim_r, pred_maskim_g, pred_maskim_b);
    
%     figure(3);
%     imshowpair(gt_maskim, pred_maskim, 'montage');
%     waitforbuttonpress;
    % save the output images
    imwrite(gt_maskim, fullfile(output_imdir, strrep(tempimname, '.tif', '_gt_overlay2.tif')));
    imwrite(pred_maskim, fullfile(output_imdir, strrep(tempimname, '.tif', '_unetpred_overlay2.tif')));
    
end
