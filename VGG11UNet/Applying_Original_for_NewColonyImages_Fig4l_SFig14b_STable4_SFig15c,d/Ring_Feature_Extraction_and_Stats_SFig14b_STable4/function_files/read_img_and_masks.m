%% read_img_and_masks

function [this_img, pred_mask, skel_mask] = read_img_and_masks(imgsFolder, predMasksFolder, skelMasksFolder, img_name)

    img_path = strcat(imgsFolder,'/',img_name);
    this_img = im2double(imread(img_path));
    
    pred_name = insertBefore(img_name, '.tif', '_pred');
    pred_mask_path = strcat(predMasksFolder,'/',pred_name);
    pred_mask = im2double(imread(pred_mask_path));
    
    skel_name = insertBefore(img_name, '.tif', '_skel');
    skel_mask_path = strcat(skelMasksFolder,'/',skel_name);
    skel_mask = im2double(imread(skel_mask_path));

end
