%% red_masks_and_overlays

function [redPred, redFin, redPredOver, redFinOver_Yel] = red_masks_and_overlays(pred_mask, cleaned_mask, this_img, final_cols_to_consider)

    % red color map
    cmap = [1 0 0];
    
    % red version of u-net prediction
    redPred = label2rgb(logical(pred_mask),cmap);

    % red version of final post-processed mask
    redFin = label2rgb(logical(cleaned_mask),cmap);
    
    % overlay red predicted mask on polarim
    redPredOver = imoverlay(this_img,pred_mask,cmap);
    % overlay red final post-processed mask on polarim
    redFinOver = imoverlay(this_img,cleaned_mask,cmap);
    
    % highlight ignored cols in yellow & make final overlay
    col_mask = logical(zeros(size(cleaned_mask)));
    all_cols = 1:1000;
    omit_cols = setdiff(all_cols,final_cols_to_consider);
    col_mask(:,omit_cols) = 1;
    ymap = [.99 .91 .51]; % yellow crayola
    redFinOver_Yel = labeloverlay(redFinOver,col_mask,...
                              'Colormap',ymap,'Transparency',0.6);


end