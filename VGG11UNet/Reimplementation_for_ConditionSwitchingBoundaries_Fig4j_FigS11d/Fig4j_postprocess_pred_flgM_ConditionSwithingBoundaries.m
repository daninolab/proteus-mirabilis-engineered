%% postprocess_pred_flgM_ConditionSwithingBoundaries
% Script for post-processing the masks predicted by the
% re-implemented/re-trained VGG-11 U-Net that delineate ring
% boundaries on flgM images where condition switched occurred 
% (i.e. between 25C benchtop & 37C incubator)

clear all;
clc;

% List colony images (polar forms) for which we also have predicted masks
% (all necessary images placed on current path)
img_names = {'07-20-21_flgM_5i_B_ben24_ben38_inc62_002_polarim.tif',...
             'flgM_i_5_bibi_75h_004_polarim.tif'};

% define the corresponding predicted mask names  
% (we used the model obtained after epoch 60 to make predictions)
% (aka epoch 59, where the 1st epoch is epoch 0)
% (predictions were first skeletonized in the original cloud-based...  
    % ...notebook for model implementation/evaluation)
epoch_name = '_skel_ep59'; 
mask_names = [];
for name_i = img_names
    mask_name_i = insertBefore(name_i,'.tif',epoch_name);
    mask_names = [mask_names mask_name_i];
end

% Define some variables for basic morphological operations.

% For dilation to ensure a given boundary is fully connected: 
se_di_line = strel('line',22,0); % for horizontal connectivity 
se_di_disk = strel('disk',2); % for vertical+horizontal connectivity 

% For subsequent skeletonizing:
branch_len = 100;

% For connecting boundaries to left & right edges of image matrix
indent = 24; % approx. max index inward from edges
se_di_edge = strel('line',indent*2,0); 

for n = 1:length(img_names)
    
    % read in colony image and corresponding pre-skeletonized, predicted mask
    i_name = img_names{n};
    img = imread(i_name);
    m_name = mask_names{n};
    mask = imread(m_name);
    
    % first visualize pre-skeletonized prediction on colony image
    orig_over = imoverlay(img,mask,'r');
    figure();
    imshow(orig_over);
    title('Pre-skeletonized, original prediction');
    
    % dilation to ensure a given boundary is fully connected
    di_mask_1 = imdilate(mask,se_di_line);
    di_over_1 = imoverlay(img,di_mask_1,'r');
    figure();
    imshow(di_over_1);
    title('Horizontally dilated prediction');
    
    di_mask_2 = imdilate(di_mask_1,se_di_disk);
    di_over_2 = imoverlay(img,di_mask_2,'r');
    figure();
    imshow(di_over_2);
    title('Spherically dilated prediction');
    pause(.2);
    
    % skeletonize to obtain single-pixel thick boundaries
    skel_mask = bwskel(logical(di_mask_2),'MinBranchLength',branch_len);
    skel_over = imoverlay(img,skel_mask,'r');
    figure();
    imshow(skel_over);
    title('Skeletonized, dilated prediction');
    
    % apply flat line-shaped structuring element to dilate near the
    % left & right image edges to re-connect the boundaries with the edges
    left = imdilate(skel_mask(:,1:indent),se_di_edge);
    right = imdilate(skel_mask(:,end-indent:end),se_di_edge);
    edge_mask = skel_mask;
    edge_mask(:,1:indent) = left;
    edge_mask(:,end-indent:end) = right;
    edge_over = imoverlay(img,edge_mask,'r');
    figure();
    imshow(edge_over);
    title('Edge-connected prediction');
    
    % now thin edges back to 1 pixel thick
    fin_mask = edge_mask;
    
    [left_rows, left_cols] = find(edge_mask(:,1:indent));
    uniq_left_rows = unique(left_rows);
    dist_left_rows = diff(uniq_left_rows);
    adj_left_ids = find(dist_left_rows==1) + 1;
    
    for i_left = adj_left_ids
        rem_row_left = uniq_left_rows(i_left);
        fin_mask(rem_row_left,1:indent) = 0;
    end
    
    [right_rows, right_cols] = find(edge_mask(:,end-indent:end));
    uniq_right_rows = unique(right_rows);
    dist_right_rows = diff(uniq_right_rows);
    adj_right_ids = find(dist_right_rows==1) + 1;
    
    for i_right = adj_right_ids
        rem_row_right = uniq_right_rows(i_right);
        fin_mask(rem_row_right,end-indent:end) = 0;
    end
    
    fin_over = imoverlay(img,fin_mask,'r');
    figure();
    imshow(fin_over);
    title('Final post-processed prediction');
    
    % save final post-processed predicted mask to current path
    fin_name = insertBefore(m_name,'.tif','_postproc');
    imwrite(fin_mask, fin_name);
    
end

