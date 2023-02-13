%% reduce_mask_for_calcs

function [reduced_mask, vertical_count] = reduce_mask_for_calcs(num_boundaries_mode, num_cols_consider, cols_to_consider, labeled_mask)

    vertical_count = 0;
    
    [h w] = size(labeled_mask);
    reduced_mask = zeros(h,w);

    for i = 1:num_cols_consider

        % grab column
        col_idx = cols_to_consider(i);
        this_col = labeled_mask(:,col_idx);

        % see if the same ring boundary crosses a given column more than once
        % if so, only keep the lowest idx for calculations
        col_labels = unique(this_col);
        col_labels(col_labels == 0) = []; % drop the 0 (background label)

        for l = 1:length(col_labels)
            label_l = col_labels(l);
            label_l_idx = find(this_col == label_l);
            label_l_counts = length(label_l_idx);

            if label_l_counts > 1
                vertical_count = vertical_count + label_l_counts;
                lowest_idx = max(label_l_idx);
                ignore_idx = label_l_idx(label_l_idx ~= lowest_idx);
                this_col(ignore_idx) = 0;
            end
        end
        
        % insert reduced column into reduced_mask
        reduced_mask(:,col_idx) = this_col;
    end


end